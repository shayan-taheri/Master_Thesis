% *** Performing Multiple Linear Regression ***
% MFO ~ Temperature + Transformer Voltage + Fitted MFO (using Time)
%% ---- Preparing Environment ----
clear all;
close all;

card = 'b5c7'; % Name of Card

tRec = 10000; % Number of Records

testDir1 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir1 = ['/media/SHAYAN_HDD/Data/collection_6/',card];
figDir = '/media/SHAYAN_HDD/Results/Analysis_7/mfo_temp_vtran_yfit/col6';

% ---- Dataset Variables ----
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
resid_dis = y - fit_dis;

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

resid_chg = y - fit_chg;

%% ---- Temperature and Voltage Fitting ----

% ---- For Discharging Trend ----
x_dat = [x,z1,fit_dis];
mdl = LinearModel.fit(x_dat,y);
ytot_dis = feval(mdl,x_dat);
tres_dis = y - ytot_dis;

% ---- For Charging Trend ----
x_dat = [x,z1,fit_chg];
mdl = LinearModel.fit(x_dat,y);
ytot_chg = feval(mdl,x_dat);
tres_chg = y - ytot_chg;

%% ---- R-Square Calculation ----

% ---- For Discharging Trend ----
SSresid = sum(resid_dis.^2);
rsq_dis = 1 - (SSresid/SStotal);

SSresid = sum(tres_dis.^2);
trsq_dis = 1 - (SSresid/SStotal);

% ---- For Charging Trend ----
SSresid = sum(resid_chg.^2);
rsq_chg = 1 - (SSresid/SStotal);

SSresid = sum(tres_chg.^2);
trsq_chg = 1 - (SSresid/SStotal);

%% ---- Plotting Results ----

% ---- For Discharging Trend ----

fig_id = figure();
plot(resid_dis,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Time Regression', 'FontSize', 14);
title(['Time Regression for Discharging Trend',' -- ','First R^2 = ',num2str(rsq_dis)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_resid_dis_tran','.pdf']);

fig_id = figure();
plot(tres_dis,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Total Regression', 'FontSize', 14);
title(['Total Regression for Discharging Trend',' -- ','Second R^2 = ',num2str(trsq_dis)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_tres_dis_tran','.pdf']);

fig_id = figure();
plot(fit_dis);
xlabel('Index of Record', 'FontSize', 14);
ylabel('Fitted MFO of Time Regression', 'FontSize', 14);
title('Time Regression Analysis for Discharging Trend', 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_fit_dis_tran','.pdf']);

fig_id = figure();
plot(ytot_dis);
xlabel('Index of Record', 'FontSize', 14);
ylabel('Fitted MFO of Total Regression', 'FontSize', 14);
title('Total Regression Analysis for Discharging Trend', 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_ytot_dis_tran','.pdf']);

% ---- For Charging Trend ----

fig_id = figure();
plot(resid_chg,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Time Regression', 'FontSize', 14);
title(['Time Regression for Charging Trend',' -- ','First R^2 = ',num2str(rsq_chg)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_resid_chg_tran','.pdf']);

fig_id = figure();
plot(tres_chg,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Total Regression', 'FontSize', 14);
title(['Total Regression for Charging Trend',' -- ','Second R^2 = ',num2str(trsq_chg)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_tres_chg_tran','.pdf']);

fig_id = figure();
plot(fit_chg);
xlabel('Index of Record', 'FontSize', 14);
ylabel('Fitted MFO of Time Regression', 'FontSize', 14);
title('Time Regression Analysis for Charging Trend', 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_fit_chg_tran','.pdf']);

fig_id = figure();
plot(ytot_chg);
xlabel('Index of Record', 'FontSize', 14);
ylabel('Fitted MFO of Total Regression', 'FontSize', 14);
title('Total Regression Analysis for Charging Trend', 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[figDir,'/',card,'_ytot_chg_tran','.pdf']);
% ----------------------------------