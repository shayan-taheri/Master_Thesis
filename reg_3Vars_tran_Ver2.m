% *** Performing Multiple Linear Regression ***
% MFO ~ Temperature + Transformer Voltage + Time Model
%% ---- Preparing Environment ----
clear all;
close all;

card = 'b5c7'; % Name of Card

tRec = 10000; % Number of Records

% ---- Dataset Variables ----
testDir1 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir1 = ['/media/SHAYAN_HDD/Data/collection_6/',card];
temp_dat1 = zeros(1,tRec);
mfo_dat1 = zeros(1,tRec);
sv_Tran1 = zeros(1,tRec);
time_dat = zeros(1,tRec);

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data of Dataset ----

disp('Getting Data for Dataset:');

% ---- Getting Temperature and Voltage Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'temp','volt_dat','ts');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'temp','volt_dat','ts');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'temp','volt_dat','ts');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'temp','volt_dat','ts');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'temp','volt_dat','ts');
    end
    
    sv_Tran1(i) = volt_dat(1);
        
    time_dat(i) = ts;
    
    if (isempty(temp) == 0)
        temp_dat1(i) = temp(1);
    else
        temp_dat1(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
temp_ind1 = find (temp_dat1~=0.0);
temp_dat1 = temp_dat1(temp_ind1);

sv_Tran1 = sv_Tran1(temp_ind1);

% ---- Getting Matched Filter Output ----
load([testDir1,'/',card,'_',card,'.mat'],'m0');
mfo_dat1 = m0.genOrig;
mfo_dat1 = mfo_dat1(temp_ind1); % Using matched filter outputs of good records.

%% ---- Calculating Time Data ----

time_dat = time_dat - time_dat(1);

time_dat = time_dat(temp_ind1);

time_dat = time_dat';

%% ---- Preparing Data ----

% ---- Preparing Data for Fitting ----
x = temp_dat1;
z1 = sv_Tran1;
y = mfo_dat1;
[x,z1,y] = prepareSurfaceData(x,z1,y);

% ---- Removing Outliers ----
[x,~] = outliers(x,10);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);
z1 = z1(~nan_x);
time_dat = time_dat(~nan_x);

[y,~] = outliers(y,10);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
z1 = z1(~nan_y);
time_dat = time_dat(~nan_y);

[z1,~] = outliers(z1,10);
nan_z1 = isnan(z1);
z1 = z1(~nan_z1);
y = y(~nan_z1);
x = x(~nan_z1);
time_dat = time_dat(~nan_z1);

SStotal = (length(y)-1) * var(y); % Total Sum of Squares

%% ---- Time Fitting ----

% ---- Discharging Trend (exp(-Time/Constant)) ----
dis_obj = fit(-time_dat,y,'exp1');
fit_dis = feval(dis_obj,-time_dat);
resid_dis1 = y - fit_dis;

% ---- Charging Trend (1 - exp(-Time/Constant)) ----
sat_beg = floor(0.75 * length(y)); % Index for Beginning of Saturation
sat_val = mean(y(sat_beg:end)); % Finding the value of saturation
tau_val = 0.63*sat_val;

for j=0:7
    y_mean(j+1) = mean(y(j*1000+1:j*1000+1000));
end

if length(y) <= 9000
    y_mean(9) = mean(y(8001:end));
else
    y_mean(9) = mean(y(8001:9000));
    y_mean(10) = mean(y(9001:end));
end

[~,tau_sec] = min(abs(y_mean-tau_val));

tau_beg = (tau_sec-1)*1000;

[~,tau_index] = min(abs(y(tau_beg+1:tau_beg+1000)-tau_val));

tau_const = time_dat(tau_index);

X_new = (1 - exp(-time_dat/tau_const));

chg_obj = polyfit(X_new,y,1);

fit_chg = polyval(chg_obj,X_new);

resid_chg1 = y - fit_chg;

%% ---- Temperature and Voltage Fitting ----

% ---- For Discharging Trend ----
Fit_Object = fit([x,z1],resid_dis1,'poly11');
p_vals = coeffvalues(Fit_Object);
yfit = p_vals(1) + p_vals(2)*x + p_vals(3)*z1;
resid_dis2 = resid_dis1 - yfit;
resid_dis3 = resid_dis1 + yfit;

% ---- For Charging Trend ----
Fit_Object = fit([x,z1],resid_chg1,'poly11');
p_vals = coeffvalues(Fit_Object);
yfit = p_vals(1) + p_vals(2)*x + p_vals(3)*z1;
resid_chg2 = resid_chg1 - yfit;
resid_chg3 = resid_chg1 + yfit;

%% ---- R-Square Calculation ----

% ---- For Discharging Trend ----
SSresid = sum(resid_dis1.^2);
rsq_dis1 = 1 - (SSresid/SStotal);

SSresid = sum(resid_dis2.^2);
rsq_dis2 = 1 - (SSresid/SStotal);

SSresid = sum(resid_dis3.^2);
rsq_dis3 = 1 - (SSresid/SStotal);

% ---- For Charging Trend ----
SSresid = sum(resid_chg1.^2);
rsq_chg1 = 1 - (SSresid/SStotal);

SSresid = sum(resid_chg2.^2);
rsq_chg2 = 1 - (SSresid/SStotal);

SSresid = sum(resid_chg3.^2);
rsq_chg3 = 1 - (SSresid/SStotal);

%% ---- Plotting Results ----

% ---- For Discharging Trend ----

fig_id = figure();
plot(resid_dis1,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual 1', 'FontSize', 14);
title(['MLR Analysis for Discharging Trend',' -- ','First R^2 = ',num2str(rsq_dis1)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_5/col6/',card,'_dis','_tran_1','.pdf']);

fig_id = figure();
plot(resid_dis2,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual 2', 'FontSize', 14);
title(['MLR Analysis for Discharging Trend',' -- ','Second R^2 = ',num2str(rsq_dis2)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_5/col6/',card,'_dis','_tran_2','.pdf']);

fig_id = figure();
plot(resid_dis3,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual 3', 'FontSize', 14);
title(['MLR Analysis for Discharging Trend',' -- ','Third R^2 = ',num2str(rsq_dis3)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_5/col6/',card,'_dis','_tran_3','.pdf']);

% ---- For Charging Trend ----

fig_id = figure();
plot(resid_chg1,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual 1', 'FontSize', 14);
title(['MLR Analysis for Charging Trend',' -- ','First R^2 = ',num2str(rsq_chg1)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_5/col6/',card,'_chg','_tran_1','.pdf']);

fig_id = figure();
plot(resid_chg2,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual 2', 'FontSize', 14);
title(['MLR Analysis for Charging Trend',' -- ','Second R^2 = ',num2str(rsq_chg2)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_5/col6/',card,'_chg','_tran_2','.pdf']);

fig_id = figure();
plot(resid_chg3,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual 3', 'FontSize', 14);
title(['MLR Analysis for Charging Trend',' -- ','Third R^2 = ',num2str(rsq_chg3)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_5/col6/',card,'_chg','_tran_3','.pdf']);
% ----------------------------------