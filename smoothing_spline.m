% *** Performing Regression using Smoothed Data ***

% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c7';
testDir = '/media/SHAYAN_HDD/Results/Collection_5/test_res_reg/1';
recsDir = ['/media/SHAYAN_HDD/Data/collection_5/',card];
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

% ---- Removing Extreme Outliers ----
%{
[x,~] = outliers(x,1);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,1);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
%}

% ---- Preparing Data ----
x = x';
y = y';
% ---- Smoothing Data ----
mdl = fit(x,y,'smoothingspline');
ysmooth = feval(mdl,x);

% ---- Linear Regression 1 ----
[sdf1,yinf1] = fit(x,ysmooth,'poly1','Robust','LAR');
yfit1 = feval(sdf1,x);
% ---- R-Square Calculation 1 ----
rsq1 = yinf1.rsquare;

% ---- Linear Regression 2 ----
[sdf2,yinf2] = fit(x,ysmooth,'poly1','Robust','Bisquare');
yfit2 = feval(sdf2,x);
% ---- R-Square Calculation 2 ----
rsq2 = yinf2.rsquare;

% ---- Plotting Results ----
fig_id = figure();

subplot(2,1,1);
plot(x,y,'*b');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title('Robust Regression on Smoothed Data','FontSize',14);
set(gca, 'fontsize', 12);

subplot(2,1,2);
plot(x,ysmooth,'pm');
hold on;
plot(x,yfit1,'-k');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title(['R^2 of Least Absolute Residuals Method = ',num2str(rsq1)], 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Smoothed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',7);
hold off;

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_2/smooth/col5/',card,'_lar','.pdf']);

fig_id = figure();

subplot(2,1,1);
plot(x,y,'*b');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title('Robust Regression on Smoothed Data','FontSize',14);
set(gca, 'fontsize', 12);

subplot(2,1,2);
plot(x,ysmooth,'pc');
hold on;
plot(x,yfit2,'-k');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title(['R^2 of Bisquare Weights Method = ',num2str(rsq2)], 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Smoothed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',7);
hold off;

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_2/smooth/col5/',card,'_bw','.pdf']);
% -------------------------------