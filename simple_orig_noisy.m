% *** Performing Linear Regression using Noisy Sample Points ***


%% ---- Preparing Environment ----

clear all;
close all;
deg = 1; % Degree of Regression
card = 'b5c7';
tSamp = 37900; % Number of Noisy Sample Points

testDir = '/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name';
figDir = '/media/SHAYAN_HDD/Results/Analysis/noisy/col2';
recsDir = ['/media/SHAYAN_HDD/Data/collection_2/',card];
tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Regression Analysis ----

% ---- Fetching Records Data ----

x1 = zeros(1,tRec); % Vector for Mean of Noisy Sample Points
x2 = zeros(1,tRec); % Vector for Standard Deviation of Noisy Sample Points

for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir,'/sample','0000',num2str(i),'.mat'],'rec');
    elseif (i >= 10) && (i <= 99)
        load([recsDir,'/sample','000',num2str(i),'.mat'],'rec');
    elseif (i >= 100) && (i <= 999)
        load([recsDir,'/sample','00',num2str(i),'.mat'],'rec');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir,'/sample','0',num2str(i),'.mat'],'rec');
    else
        load([recsDir,'/sample',num2str(i),'.mat'],'rec');
    end
    
    x1(i) = mean(rec(1:tSamp));
    
    x2(i) = std(rec(1:tSamp));
    
end

% ---- Fetching Matched Filter Output ----
y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;

%% ---- Performing Linear Regression on Mean Data ----

% ---- Linear Regression ----
p = polyfit(x1,y,deg);
yfit = polyval(p,x1);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p));
% ---- Plotting Results ----
fig_id = figure();
plot(x1,y,'*b');
hold on;
plot(x1,yfit,'-r');
xlabel('Mean of Noisy Sample Points', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','NorthEast');
set(fig_leg,'FontSize',10);
hold off;
% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);
% ---- Saving Results ----
saveas(fig_id,[figDir,'/',card,'_','mean_',num2str(deg),'.pdf']);

%% ---- Performing Linear Regression on Standard Deviation Data ----

% ---- Linear Regression ----
p = polyfit(x2,y,deg);
yfit = polyval(p,x2);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p));
% ---- Plotting Results ----
fig_id = figure();
plot(x2,y,'*b');
hold on;
plot(x2,yfit,'-r');
xlabel('Standard Deviation of Noisy Sample Points', 'FontSize', 12);
ylabel('Matched Filter Output', 'FontSize', 12);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','NorthEast');
set(fig_leg,'FontSize',10);
hold off;
% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);
% ---- Saving Results ----
saveas(fig_id,[figDir,'/',card,'_','std_',num2str(deg),'.pdf']);

% ----------------------------------------------------