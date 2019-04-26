% *** MFO ~ Temp. + Supply Voltage + Mean/Standard Deviation of Noise ***
%% ---- Preparing Environment ----
clear all;
close all;

card = 'b5c7'; % Name of Card

tRec = 10000; % Number of Records

tSamp = 37900; % Number of Noisy Sample Points

% ---- Dataset 1 Variables ----
testDir1 = '/media/SHAYAN_HDD/Results/Collection_5/test_res_reg/1';
recsDir1 = ['/media/SHAYAN_HDD/Data/collection_5/',card];
temp_dat1 = zeros(1,tRec);
mfo_dat1 = zeros(1,tRec);
sv_Tran1 = zeros(1,tRec);
sv_IC1 = zeros(1,tRec);
noisy_mean1 = zeros(1,tRec);
noisy_std1 = zeros(1,tRec);

% ---- Dataset 2 Variables ----
testDir2 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir2 = ['/media/SHAYAN_HDD/Data/collection_6/',card];
temp_dat2 = zeros(1,tRec);
mfo_dat2 = zeros(1,tRec);
sv_Tran2 = zeros(1,tRec);
sv_IC2 = zeros(1,tRec);
noisy_mean2 = zeros(1,tRec);
noisy_std2 = zeros(1,tRec);

% ---- Directory of Results ----
figDir = '/media/SHAYAN_HDD/Results/Analysis_6/mfo_temp_volt_noise';

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data of Dataset 1 ----

disp('Getting Data for Dataset 1:');

% ---- Getting Temperature and Voltage Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'temp','volt_dat','rec');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'temp','volt_dat','rec');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'temp','volt_dat','rec');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'temp','volt_dat','rec');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'temp','volt_dat','rec');
    end
    
    if (isempty(volt_dat) == 0)
        sv_Tran1(i) = volt_dat(1);
        sv_IC1(i) = volt_dat(2);
    else
        sv_Tran1(i) = 0.0;
        sv_IC1(i) = 0.0;
    end
    
    if (isempty(temp) == 0)
        temp_dat1(i) = temp(1);
    else
        temp_dat1(i) = 0.0;
    end
    
    noisy_mean1(i) = mean(rec(1:tSamp));
    
    noisy_std1(i) = std(rec(1:tSamp));
    
end

% Finding Good Records: The ones that have data.
temp_ind1 = find (temp_dat1~=0.0);
temp_dat1 = temp_dat1(temp_ind1);

sv_Tran1 = sv_Tran1(temp_ind1);

sv_IC1 = sv_IC1(temp_ind1);

noisy_mean1 = noisy_mean1(temp_ind1);
    
noisy_std1 = noisy_std1(temp_ind1);

% ---- Getting Matched Filter Output ----
load([testDir1,'/',card,'_',card,'.mat'],'m0');
mfo_dat1 = m0.genOrig;
mfo_dat1 = mfo_dat1(temp_ind1); % Using matched filter outputs of good records.

%% ---- Getting Data of Dataset 2 ----

disp('Getting Data for Dataset 2:');

% ---- Getting Temperature and Voltage Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir2,'/sample','0000',num2str(i),'.mat'],'temp','volt_dat','rec');
    elseif (i >= 10) && (i <= 99)
        load([recsDir2,'/sample','000',num2str(i),'.mat'],'temp','volt_dat','rec');
    elseif (i >= 100) && (i <= 999)
        load([recsDir2,'/sample','00',num2str(i),'.mat'],'temp','volt_dat','rec');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir2,'/sample','0',num2str(i),'.mat'],'temp','volt_dat','rec');
    else
        load([recsDir2,'/sample',num2str(i),'.mat'],'temp','volt_dat','rec');
    end
    
    if (isempty(volt_dat) == 0)
        sv_Tran2(i) = volt_dat(1);
        sv_IC2(i) = volt_dat(2);
    else
        sv_Tran2(i) = 0.0;
        sv_IC2(i) = 0.0;
    end
    
    if (isempty(temp) == 0)
        temp_dat2(i) = temp(1);
    else
        temp_dat2(i) = 0.0;
    end
    
    noisy_mean2(i) = mean(rec(1:tSamp));
    
    noisy_std2(i) = std(rec(1:tSamp));
    
end

% Finding Good Records: The ones that have data.
temp_ind2 = find (temp_dat2~=0.0);
temp_dat2 = temp_dat2(temp_ind2);

sv_Tran2 = sv_Tran2(temp_ind2);

sv_IC2 = sv_IC2(temp_ind2);

noisy_mean2 = noisy_mean2(temp_ind2);
    
noisy_std2 = noisy_std2(temp_ind2);

% ---- Getting Matched Filter Output ----
load([testDir2,'/',card,'_',card,'.mat'],'m0');
mfo_dat2 = m0.genOrig;
mfo_dat2 = mfo_dat2(temp_ind2); % Using matched filter outputs of good records.

%% ---- Preparing Data and Removing Outliers ----

% ---- Preparing Data for Fitting ----
mfo_tot = [mfo_dat1,mfo_dat2]';
temp_tot = [temp_dat1,temp_dat2]';
tran_tot = [sv_Tran1,sv_Tran2]';
ic_tot = [sv_IC1,sv_IC2]';
mean_tot = [noisy_mean1,noisy_mean2]';
std_tot = [noisy_std1,noisy_std2]';

% ---- Removing Outliers ----
[mfo_tot,~] = outliers(mfo_tot,5);
nan_ind = isnan(mfo_tot);

mfo_tot = mfo_tot(~nan_ind);
temp_tot = temp_tot(~nan_ind);
tran_tot = tran_tot(~nan_ind);
mean_tot = mean_tot(~nan_ind);
std_tot = std_tot(~nan_ind);
ic_tot = ic_tot(~nan_ind);
% ---------------------------
[temp_tot,~] = outliers(temp_tot,5);
nan_ind = isnan(temp_tot);

mfo_tot = mfo_tot(~nan_ind);
temp_tot = temp_tot(~nan_ind);
tran_tot = tran_tot(~nan_ind);
mean_tot = mean_tot(~nan_ind);
std_tot = std_tot(~nan_ind);
ic_tot = ic_tot(~nan_ind);
% ---------------------------
[mean_tot,~] = outliers(mean_tot,5);
nan_ind = isnan(mean_tot);

mfo_tot = mfo_tot(~nan_ind);
temp_tot = temp_tot(~nan_ind);
tran_tot = tran_tot(~nan_ind);
mean_tot = mean_tot(~nan_ind);
std_tot = std_tot(~nan_ind);
ic_tot = ic_tot(~nan_ind);
% ---------------------------
[std_tot,~] = outliers(std_tot,5);
nan_ind = isnan(std_tot);

mfo_tot = mfo_tot(~nan_ind);
temp_tot = temp_tot(~nan_ind);
tran_tot = tran_tot(~nan_ind);
mean_tot = mean_tot(~nan_ind);
std_tot = std_tot(~nan_ind);
ic_tot = ic_tot(~nan_ind);

% ---- Storing Major Data ----
mfo_rec = mfo_tot;
temp_rec = temp_tot;
mean_rec = mean_tot;
std_rec = std_tot;

%% ---- Regression Analysis using Transformer Supply Voltage ----

% ---- Removing Outliers of Transformer Supply Voltage ----
[tran_tot,~] = outliers(tran_tot,5);
nan_ind = isnan(tran_tot);

mfo_tot = mfo_tot(~nan_ind);
temp_tot = temp_tot(~nan_ind);
tran_tot = tran_tot(~nan_ind);
mean_tot = mean_tot(~nan_ind);
std_tot = std_tot(~nan_ind);

% ---- Performing Multiple Linear Regression 1 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,tran_tot,mean_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y);
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 1',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr1','_tran','.pdf']);

% ---- Performing Multiple Linear Regression 2 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,tran_tot,std_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y);
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 2',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr2','_tran','.pdf']);

% ---- Performing Multiple Linear Regression 3 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,tran_tot,mean_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 3',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr3','_tran','.pdf']);

% ---- Performing Multiple Linear Regression 4 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,tran_tot,std_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 4',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr4','_tran','.pdf']);

% ---- Performing Multiple Linear Regression 5 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,tran_tot,(std_tot ./ mean_tot)];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 5',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr5','_tran','.pdf']);

% ---- Performing Multiple Linear Regression 6 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,tran_tot,(std_tot .* mean_tot)];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 6',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr6','_tran','.pdf']);

%% ---- Regression Analysis using IC Supply Voltage ----

% ---- Restoring Major Data ----
mfo_tot = mfo_rec;
temp_tot = temp_rec;
mean_tot = mean_rec;
std_tot = std_rec;

% ---- Removing Outliers of IC Supply Voltage ----
[ic_tot,~] = outliers(ic_tot,5);
nan_ind = isnan(ic_tot);

mfo_tot = mfo_tot(~nan_ind);
temp_tot = temp_tot(~nan_ind);
ic_tot = ic_tot(~nan_ind);
mean_tot = mean_tot(~nan_ind);
std_tot = std_tot(~nan_ind);

% ---- Performing Multiple Linear Regression 1 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,ic_tot,mean_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y);
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 1',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr1','_ic','.pdf']);

% ---- Performing Multiple Linear Regression 2 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,ic_tot,std_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y);
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 2',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr2','_ic','.pdf']);

% ---- Performing Multiple Linear Regression 3 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,ic_tot,mean_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 3',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr3','_ic','.pdf']);

% ---- Performing Multiple Linear Regression 4 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,ic_tot,std_tot];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 4',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr4','_ic','.pdf']);

% ---- Performing Multiple Linear Regression 5 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,ic_tot,(std_tot ./ mean_tot)];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 5',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr5','_ic','.pdf']);

% ---- Performing Multiple Linear Regression 6 ----

% ---- Estimating MFO ----
x_dat = [temp_tot,ic_tot,(std_tot .* mean_tot)];
y = mfo_tot;
mdl = LinearModel.fit(x_dat,y,'RobustOpts','on');
yfit = feval(mdl,x_dat);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title(['MLR Analysis 6',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_mlr6','_ic','.pdf']);
% ----------------------------------------------------