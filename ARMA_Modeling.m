% *** Autoregressive Moving Average (ARMA) Modeling ***

% ---- Wikipedia Equation ----
% Xt = e,t + Sum(phi,i * X,t-i) + Sum(theta,i * e,t-i)

% ---- MATLAB Equation ----
% A(q) * y(t) = C(q) * e(t)

%% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c7';

testDir = '/media/SHAYAN_HDD/Results/Collection_7/test_res_reg/1'; % Directory of MFO

resDir = '/media/SHAYAN_HDD/Results/Analysis_12'; % Directory of Results

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

% ---- Matched Filter Output ----
mfo = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo = m0.genOrig;

%% ---- Preparing Variables ----
y = mfo'; % Transpose MFO Data Due to ARMAX Function.

% ---- Removing Outliers ----
[y,~] = outliers(y,30);
nan_y = isnan(y);
y = y(~nan_y);

%% ---- Performing ARMA Modeling ----

% ---- Data Object Creation ----
data = iddata(y,[],1);

% ---- Best Value for Goodness of Fit ----
best_fit = 0.0;
best_AR = 0.0;
best_MA = 0.0;

% ---- The Values for Orders ----
Order_Values = [1:5,10,20,50,100,110]; % 10 Numbers --> 10^2 Possibilities

% ---- Performing ARMA Modeling ----

disp('Calculating Results:');

for j1 = 1:size(Order_Values,2)

	AR_Order = Order_Values(j1); % Order of Autoregression
    
	for j2 = 1:size(Order_Values,2)
        
        MA_Order = Order_Values(j2); % Order of Moving Average
        
        try
            
            system = armax(data,[AR_Order,MA_Order]);

            [~,fit_value,~] = compare(data,system);
            
        catch
            
            continue;
            
        end

        if (best_fit < fit_value) % Getting the Best Fit.

            best_fit = fit_value;
            best_AR = AR_Order;
            best_MA = MA_Order;

        end

	end
    
	if mod(j1,1) == 0
        
        disp(['Percentage of Completion = ',num2str((j1 * 100)/size(Order_Values,2)),'%']);
        
	end

end

% ---- Plotting Results ----
system = armax(data,[best_AR,best_MA]);
[~,fit_value,~] = compare(data,system);

fig_id = figure();
compare(data,system);
saveas(fig_id,[resDir,'/',card,'_','arma_fit','.pdf']);

fig_id = figure();
resid(system,data);
xlabel('Lag', 'FontSize', 14);
ylabel('Autocorrelation Output', 'FontSize', 14);
title(['R-Square Value (%) = ',num2str(fit_value),' , AR Order = ',num2str(best_AR), ...
    ' , MA Order = ',num2str(best_MA)], 'FontSize', 12);
set(gca, 'fontsize', 12);

saveas(fig_id,[resDir,'/',card,'_','arma_res','.pdf']);
% -------------------------------