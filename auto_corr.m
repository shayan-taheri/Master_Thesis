% *** Performing Auto Correlation using Original Data ***
% ---- Preparing Environment ----
clear all;
close all;
deg = 1; % Degree of Fitting Data
card = 'b5c7';
recsDir = '/media/SHAYAN_HDD/Data/collection_2/b5c7';
testDir = '/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name';
tRec = 10000; % Number of Records
% ---- Fetching Temperature Data ----
x = zeros(1,tRec); % Vector for Temperature Data
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir,'/sample',num2str(i),'.mat'],'temp');
    end
    
    if (isempty(temp) == 0)
        x(i) = temp(1);
    else
        x(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
g_ix = find (x~=0.0);
x = x(g_ix);

% ---- Matched Filter Output ----
y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;
y = y(g_ix); % Using matched filter outputs of good records.

% ---- Removing Outliers ----
[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

% ---- Residual 1 ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);
yresid1 = y - yfit;
% ---- Residual 2 ----
y_hat = y/max(y);
x_hat = x/max(x);
yresid2 = y_hat - x_hat;
% ---- Plotting Results ----
fig_id = figure();
autocorr(yresid1);
xlabel('Lag', 'FontSize', 12);
ylabel('Sample Autocorrelation', 'FontSize', 12);
title('Residual = Actual MFO - Predicted MFO', 'FontSize', 12);
set(gca, 'fontsize', 12);
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_2/auto_corr/col2/',card,'_','res1','.pdf']);

fig_id = figure();
autocorr(yresid2);
xlabel('Lag', 'FontSize', 12);
ylabel('Sample Autocorrelation', 'FontSize', 12);
title('Residual = Relative MFO - Relative Temperature', 'FontSize', 12);
set(gca, 'fontsize', 12);
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_2/auto_corr/col2/',card,'_','res2','.pdf']);
% -------------------------------