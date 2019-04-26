% *** Performing Multiple Linear Regression ***
% MFO ~ Temperature + IC Voltage + (1 - exp(-Time/Constant))
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
sv_IC1 = zeros(1,tRec);
date_raw = zeros(1,tRec);
date_diff = zeros(1,tRec);

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
    
    sv_IC1(i) = volt_dat(2);
        
    date_raw(i) = ts;
    
    if (isempty(temp) == 0)
        temp_dat1(i) = temp(1);
    else
        temp_dat1(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
temp_ind1 = find (temp_dat1~=0.0);
temp_dat1 = temp_dat1(temp_ind1);

sv_IC1 = sv_IC1(temp_ind1);

% ---- Getting Matched Filter Output ----
load([testDir1,'/',card,'_',card,'.mat'],'m0');
mfo_dat1 = m0.genOrig;
mfo_dat1 = mfo_dat1(temp_ind1); % Using matched filter outputs of good records.

%% ---- Calculating Time Data ----

for j=2:length(date_raw)
    date_diff(j) = date_raw(j)-date_raw(j-1);
end

date_diff = datevec(date_diff);

date_diff(:,6) = date_diff(:,4) * 3600 + date_diff(:,5) * 60 + date_diff(:,6);

date_diff = date_diff(:,6);

date_diff(1) = mean(date_diff(2:end));

date_diff = [0.0,date_diff'];

time_dat = date_diff;

for k=2:length(time_dat)
    time_dat(k) = time_dat(k) + time_dat(k-1);
end

time_dat = time_dat(2:end);

time_dat = time_dat(temp_ind1);

%% ---- Preparing Data and Getting Residual ----

% ---- Preparing Data for Fitting ----
x = temp_dat1;
z1 = sv_IC1;
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

% ---- Fitting Prepared Data ----
Fit_Object = fit([x,z1],y,'poly11');
% ---- Getting Values of Coefficients ----
p_vals = coeffvalues(Fit_Object);
% ---- Calculating Estimated Response Data ----
yfit = p_vals(1) + p_vals(2)*x + p_vals(3)*z1;
% ---- Residual Calculation ----
yresid = y - yfit;

%% ---- Getting Time Constant for Charging and Discharging Trends ----

% ---- Discharging Trend (exp(-Time/Constant)) ----
dis_obj = fit(-time_dat',yresid,'exp1');
dis_cof = coeffvalues(dis_obj);
dis_cof = 1/dis_cof(2);

% ---- Charging Trend (1 - exp(-Time/Constant)) ----
sat_beg = floor(2/3 * length(yresid)); % Index for Beginning of Saturation
sat_val = mean(yresid(sat_beg:end)); % Finding the value of saturation
Yresid = log(1-(yresid/sat_val));
chg_cof = regress(Yresid,-time_dat');
chg_cof = 1/chg_cof;

%% ---- Performing Multiple Linear Regression (Discharge) ----

% ---- Estimating MFO ----
z2 = exp(-time_dat'/dis_cof);
x_dat = [x,z1,z2];
x_dat = real(x_dat);
mdl = LinearModel.fit(x_dat,y);
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
stem(yresid);
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis for Discharging Trend',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/col6/multiple_reg/',card,'_dis','_ic','.pdf']);

%% ---- Performing Multiple Linear Regression (Charge) ----

% ---- Estimating MFO ----
z2 = 1 - exp(-time_dat'/chg_cof);
x_dat = [x,z1,z2];
x_dat = real(x_dat);
mdl = LinearModel.fit(x_dat,y);
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
stem(yresid);
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis for Charging Trend',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/col6/multiple_reg/',card,'_chg','_ic','.pdf']);
% -------------------------------