% *** Autoregressive Moving Average with Exogenous Inputs (ARMAX) Modeling ***

% ---- Wikipedia Equation ----
% Xt = e,t + Sum(phi,i * X,t-i) + Sum(theta,i * e,t-i) + Sum(eta,i * d,t-i)

% ---- MATLAB Equation ----
% A(q) * y(t) = Sum(B,i(q) * u,i(t - (n * k,i))) + C(q) * e(t)

%% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c3';

testDir = '/media/SHAYAN_HDD/Results/Collection_7/test_res_reg/1'; % Directory of MFO

recsDir1 = ['/media/SHAYAN_HDD_2/Data/2015-SPRING-07/',card]; % Directory of Records

resDir = '/media/SHAYAN_HDD/Results/Analysis_12'; % Directory of Results

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

disp('Fetching Data for Dataset:');

% ---- Getting All Data ----

temp_vect = zeros(1,tRec); % Vector for Temperature Data

for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'temp');
    end
    
    if (isempty(temp) == 0)
        temp_vect(i) = temp(1);
    else
        temp_vect(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
temp_ind = find (temp_vect~=0.0);
temp_vect = temp_vect(temp_ind);

% ---- Matched Filter Output ----
mfo = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo = m0.genOrig;
mfo = mfo(temp_ind); % Using matched filter outputs of good records.

%% ---- Preparing Variables ----
x = temp_vect'; % Transpose Temperature Data Due to ARMAX Function.
y = mfo'; % Transpose MFO Data Due to ARMAX Function.

% ---- Removing Outliers ----
[x,~] = outliers(x,30);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,30);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

%% ---- Performing ARMAX Modeling ----

% ---- Data Object Creation ----
data = iddata(y,x,1);

% ---- Best Value for Goodness of Fit ----
best_fit = 0.0;
best_AR = 0.0;
best_X = 0.0;
best_MA = 0.0;
best_Delay = 0.0; % There is no delay between input and output (Assumption).

% ---- The Values for Orders ----
Order_Values = [1:5,10,20,50,100,110]; % 10 Numbers --> 10^3 Possibilities

% ---- Performing ARMAX Modeling ----

disp('Calculating Results:');

for j1 = 1:size(Order_Values,2)
    
    AR_Order = Order_Values(j1); % Order of Autoregression
   
    for j2 = 1:size(Order_Values,2)
        
        X_Order = Order_Values(j2); % Order of Exogenous Inputs
        
        for j3 = 1:size(Order_Values,2)
            
            MA_Order = Order_Values(j3); % Order of Moving Average
            
            try

                system = armax(data,[AR_Order,X_Order,MA_Order,best_Delay]);

                [~,fit_value,~] = compare(data,system);

            catch

                continue;

            end

            if (best_fit < fit_value) % Getting the Best Fit.

                best_fit = fit_value;
                best_AR = AR_Order;
                best_X = X_Order;
                best_MA = MA_Order;
                
            end

        end
        
    end
    
	if mod(j1,1) == 0
        
        disp(['Percentage of Completion = ',num2str((j1 * 100)/size(Order_Values,2)),'%']);
        
	end
    
end

% ---- Plotting Results ----
system = armax(data,[best_AR,best_X,best_MA,best_Delay]);
[~,fit_value,~] = compare(data,system);

fig_id = figure();
compare(data,system);
saveas(fig_id,[resDir,'/',card,'_','armax_fit','.pdf']);

fig_id = figure();
resid(system,data);
subplot(2,1,1);
xlabel('Lag', 'FontSize', 14);
ylabel('Autocorrelation Output', 'FontSize', 14);
title(['R-Square Value (%) = ',num2str(fit_value)], 'FontSize', 14);
set(gca, 'fontsize', 12);

subplot(2,1,2);
xlabel('Lag', 'FontSize', 14);
ylabel('Cross Correlation Output', 'FontSize', 14);
title(['AR Order = ',num2str(best_AR),' , MA Order = ',num2str(best_MA), ...
    ' , X Order = ',num2str(best_X),' , Delay Order = ',num2str(best_Delay)], 'FontSize', 12);
set(gca, 'fontsize', 12);

saveas(fig_id,[resDir,'/',card,'_','armax_res','.pdf']);
% -------------------------------