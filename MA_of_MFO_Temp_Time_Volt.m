% *** Implementation of the Following Equations ***

% 1) MA(MFO) ~ MA(Temperature) + MA(Time Model)

% 2) MA(MFO) ~ MA(Temperature) + MA(IC Voltage) + MA(Time Model)

% Notice 1: MA() = Moving Average --> Trend of Signal

% Notice 2: Time Model --> Charging or Discharging Behavior

%% ---- Preparing Environment ----

clear all;
close all;

wind_size = 25; % Window Size for Moving Average

card = 'b5c7';

testDir = '/media/SHAYAN_HDD/Results/Collection_8/test_res_reg/1';
recsDir1 = ['/media/SHAYAN_HDD/Data/2015-SPRING-08/',card];

resDir = '/media/SHAYAN_HDD/Results/Analysis_14/MA_of_MFO_Temp_Time_Volt/col8';

tRec = 10000; % Number of Records

% ---- Variables for Analysis ----
temp_vect = zeros(1,tRec); % Vector for Temperature Data
mfo_vect = zeros(1,tRec); % Vector for MFO Data
IC_voltage = zeros(1,tRec); % Vector for IC Supply Voltage
time_vect = zeros(1,tRec); % Vector for Time Data

%% ---- Getting Data of Dataset ----

disp('Getting Data for Dataset:');

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

% ---- Getting Temperature, Voltage, and Time Data ----
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
    
    IC_voltage(i) = volt_dat(2);
        
    time_vect(i) = ts;
    
    if (isempty(temp) == 0)
        temp_vect(i) = temp(1);
    else
        temp_vect(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
temp_ind = find(temp_vect~=0.0);
temp_vect = temp_vect(temp_ind);

IC_voltage = IC_voltage(temp_ind);

% ---- Getting Matched Filter Output ----
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo_vect = m0.genOrig;
mfo_vect = mfo_vect(temp_ind); % Using matched filter outputs of good records.

% ---- Calculating Time Data ----

time_vect = time_vect - time_vect(1);

time_vect = time_vect(temp_ind);

%% ---- Removing Outliers ----

[temp_vect,~] = outliers(temp_vect,20);
nan_ind = isnan(temp_vect);

temp_vect = temp_vect(~nan_ind);
IC_voltage = IC_voltage(~nan_ind);
mfo_vect = mfo_vect(~nan_ind);
time_vect = time_vect(~nan_ind);
% ---------------------------
[IC_voltage,~] = outliers(IC_voltage,20);
nan_ind = isnan(IC_voltage);

temp_vect = temp_vect(~nan_ind);
IC_voltage = IC_voltage(~nan_ind);
mfo_vect = mfo_vect(~nan_ind);
time_vect = time_vect(~nan_ind);
% ---------------------------
[mfo_vect,~] = outliers(mfo_vect,20);
nan_ind = isnan(mfo_vect);

temp_vect = temp_vect(~nan_ind);
IC_voltage = IC_voltage(~nan_ind);
mfo_vect = mfo_vect(~nan_ind);
time_vect = time_vect(~nan_ind);
% ---------------------------
[time_vect,~] = outliers(time_vect,20);
nan_ind = isnan(time_vect);

temp_vect = temp_vect(~nan_ind);
IC_voltage = IC_voltage(~nan_ind);
mfo_vect = mfo_vect(~nan_ind);
time_vect = time_vect(~nan_ind);

%% ---- Time Fitting ----

% ---- Making All Data in Column Order (For Fitting Data) ----
temp_vect = temp_vect';
IC_voltage = IC_voltage';
mfo_vect = mfo_vect';
time_vect = time_vect';

% ---- Discharging Trend (exp(-Time/Constant)) ----
dis_obj = fit(-time_vect,mfo_vect,'exp1');
fit_dis = feval(dis_obj,-time_vect);

% ---- Charging Trend (1 - exp(-Time/Constant)) ----
sat_beg = floor(0.75 * length(mfo_vect)); % Index for Beginning of Saturation
sat_val = mean(mfo_vect(sat_beg:end)); % Finding the value of saturation
tau_val = 0.63*sat_val;

for j=0:7
    mfo_mean(j+1) = mean(mfo_vect(j*1000+1:j*1000+1000));
end

if length(mfo_vect) <= 9000
    mfo_mean(9) = mean(mfo_vect(8001:end));
else
    mfo_mean(9) = mean(mfo_vect(8001:9000));
    mfo_mean(10) = mean(mfo_vect(9001:end));
end

[~,tau_sec] = min(abs(mfo_mean-tau_val));

tau_beg = (tau_sec-1)*1000;

if (tau_beg+1000) <= length(mfo_vect)
    
    [~,tau_index] = min(abs(mfo_vect(tau_beg+1:tau_beg+1000)-tau_val));
    
else
    
    [~,tau_index] = min(abs(mfo_vect(tau_beg+1:end)-tau_val));
    
end

tau_const = time_vect(tau_index);

X_new = (1 - exp(-time_vect/tau_const));

chg_obj = polyfit(X_new,mfo_vect,1);

fit_chg = polyval(chg_obj,X_new);

%% ---- Calculating Moving Average of Data ----

% ---- Making All Data in Row Order (For Moving Average) ----
temp_vect = temp_vect';
IC_voltage = IC_voltage';
mfo_vect = mfo_vect';
time_vect = time_vect';
fit_dis = fit_dis';
fit_chg = fit_chg';

% ---- Performing Moving Average on Data ----
temp_vect = tsmovavg(temp_vect,'s',wind_size);
nan_ind = isnan(temp_vect);
temp_vect = temp_vect(~nan_ind);
% ---------------------------
IC_voltage = tsmovavg(IC_voltage,'s',wind_size);
nan_ind = isnan(IC_voltage);
IC_voltage = IC_voltage(~nan_ind);
% ---------------------------
mfo_vect = tsmovavg(mfo_vect,'s',wind_size);
nan_ind = isnan(mfo_vect);
mfo_vect = mfo_vect(~nan_ind);
% ---------------------------
fit_dis = tsmovavg(fit_dis,'s',wind_size);
nan_ind = isnan(fit_dis);
fit_dis = fit_dis(~nan_ind);
% ---------------------------
fit_chg = tsmovavg(fit_chg,'s',wind_size);
nan_ind = isnan(fit_chg);
fit_chg = fit_chg(~nan_ind);

%% ---- Performing Multiple Linear Regression 1 ----

% ---- Making All Data in Column Order (For Fitting Data) ----
temp_vect = temp_vect';
IC_voltage = IC_voltage';
mfo_vect = mfo_vect';
time_vect = time_vect';
fit_dis = fit_dis';
fit_chg = fit_chg';

% ---- For Discharging Trend ----
Fit_Object = fit([temp_vect,fit_dis],mfo_vect,'poly11');
p_vals = coeffvalues(Fit_Object);
MLR1_dis = p_vals(1) + p_vals(2)*temp_vect + p_vals(3)*fit_dis;

% ---- For Charging Trend ----
Fit_Object = fit([temp_vect,fit_chg],mfo_vect,'poly11');
p_vals = coeffvalues(Fit_Object);
MLR1_chg = p_vals(1) + p_vals(2)*temp_vect + p_vals(3)*fit_chg;

%% ---- R-Square Calculation for MLR 1 ----

% ---- For Discharging Trend ----
yresid1_dis = mfo_vect - MLR1_dis;
SSresid = sum(yresid1_dis.^2);
SStotal = (length(mfo_vect)-1) * var(mfo_vect);
rsq1_dis = 1 - (SSresid/SStotal);

% ---- For Charging Trend ----
yresid1_chg = mfo_vect - MLR1_chg;
SSresid = sum(yresid1_chg.^2);
SStotal = (length(mfo_vect)-1) * var(mfo_vect);
rsq1_chg = 1 - (SSresid/SStotal);

%% ---- Plotting Results for MLR 1 ----

% ---- For Discharging Trend ----
fig_id = figure();
plot(yresid1_dis,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR 1 Analysis for Discharging Trend',' -- ','R^2 = ',num2str(rsq1_dis)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[resDir,'/',card,'_MLR1_dis.pdf']);

% ---- For Charging Trend ----
fig_id = figure();
plot(yresid1_chg,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR 1 Analysis for Charging Trend',' -- ','R^2 = ',num2str(rsq1_chg)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[resDir,'/',card,'_MLR1_chg.pdf']);

%% ---- Performing Multiple Linear Regression 2 ----

% ---- For Discharging Trend ----
x_dat = [temp_vect,fit_dis,IC_voltage];
y = mfo_vect;
mdl = LinearModel.fit(x_dat,y);
MLR2_dis = feval(mdl,x_dat);

% ---- For Charging Trend ----
x_dat = [temp_vect,fit_chg,IC_voltage];
y = mfo_vect;
mdl = LinearModel.fit(x_dat,y);
MLR2_chg = feval(mdl,x_dat);

%% ---- R-Square Calculation for MLR 2 ----

% ---- For Discharging Trend ----
yresid2_dis = mfo_vect - MLR2_dis;
SSresid = sum(yresid2_dis.^2);
SStotal = (length(mfo_vect)-1) * var(mfo_vect);
rsq2_dis = 1 - (SSresid/SStotal);

% ---- For Charging Trend ----
yresid2_chg = mfo_vect - MLR2_chg;
SSresid = sum(yresid2_chg.^2);
SStotal = (length(mfo_vect)-1) * var(mfo_vect);
rsq2_chg = 1 - (SSresid/SStotal);

%% ---- Plotting Results for MLR 2 ----

% ---- For Discharging Trend ----
fig_id = figure();
plot(yresid2_dis,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR 2 Analysis for Discharging Trend',' -- ','R^2 = ',num2str(rsq2_dis)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[resDir,'/',card,'_MLR2_dis.pdf']);

% ---- For Charging Trend ----
fig_id = figure();
plot(yresid2_chg,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR 2 Analysis for Charging Trend',' -- ','R^2 = ',num2str(rsq2_chg)], 'FontSize', 14);
set(gca, 'fontsize', 12);

saveas(fig_id,[resDir,'/',card,'_MLR2_chg.pdf']);
% -----------------------------------------------