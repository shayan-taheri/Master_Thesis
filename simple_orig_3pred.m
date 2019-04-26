% *** Performing Simple Linear Regression (Temperature , Supply Voltages) ***

%% ---- Preparing Environment ----
clear all;
close all;
deg = 1; % Degree of Regression
card = 'b5c7';

testDir = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';

recsDir1 = ['/media/SHAYAN_HDD/Data/collection_6/',card];

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

disp('Fetching Data for Dataset:');

% ---- Getting All Data ----

sv_Tran = zeros(1,tRec); % Vector for Supply Voltage of Transformer
sv_IC = zeros(1,tRec); % Vector for Supply Voltage of IC
temp_vect = zeros(1,tRec); % Vector for Temperature Data

for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'volt_dat','temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'volt_dat','temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'volt_dat','temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'volt_dat','temp');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'volt_dat','temp');
    end
    
    if (isempty(volt_dat) == 0)
        sv_Tran(i) = volt_dat(1);
        sv_IC(i) = volt_dat(2);
    else
        sv_Tran(i) = 0.0;
        sv_IC(i) = 0.0;
    end
    
    if (isempty(temp) == 0)
        temp_vect(i) = temp(1);
    else
        temp_vect(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
tran_ind = find (sv_Tran~=0.0);
sv_Tran = sv_Tran(tran_ind);

ic_ind = find (sv_IC~=0.0);
sv_IC = sv_IC(ic_ind);

temp_ind = find (temp_vect~=0.0);
temp_vect = temp_vect(temp_ind);

% ---- Matched Filter Output ----
mfo = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo = m0.genOrig;

%% ---- Regression Analysis for Supply Voltage of Transformer ----

% ---- Preparing Variables ----
x = sv_Tran;
y = mfo(tran_ind); % Using matched filter outputs of good records.

% ---- Removing Outliers ----
[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

% ---- Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p));

% ---- Plotting Results ----
fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Transformer Supply Voltage (V)', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_6/reg_res/',card,'_','Tran_',num2str(deg),'.pdf']);

%% ---- Regression Analysis for Supply Voltage of IC ----

% ---- Preparing Variables ----
x = sv_IC;
y = mfo(ic_ind); % Using matched filter outputs of good records.

% ---- Removing Outliers ----
[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

% ---- Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p));

% ---- Plotting Results ----
fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('IC Supply Voltage (V)', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_6/reg_res/',card,'_','IC_',num2str(deg),'.pdf']);

%% ---- Regression Analysis for Temperature ----

% ---- Preparing Variables ----
x = temp_vect;
y = mfo(temp_ind); % Using matched filter outputs of good records.

% ---- Removing Outliers ----
[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

% ---- Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p));

% ---- Plotting Results ----
fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_6/reg_res/',card,'_','temp_',num2str(deg),'.pdf']);
% -------------------------------