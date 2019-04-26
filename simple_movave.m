% *** Performing Simple Linear Regression on Moving Average of Data ***

%% ---- Preparing Environment ----

clear all;
close all;

deg = 3; % Degree of Regression
wind_size = 100; % Window Size for Moving Average

card = 'b5c7';

recsDir = '/media/SHAYAN_HDD/Data/collection_2/b5c7';
testDir = '/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name';
figDir = '/media/SHAYAN_HDD/Results/Analysis/simp_mov/res2';

tRec = 10000; % Number of Records

%% ---- Fetching Temperature Data and Getting Its Moving Average ----

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

% Calculating Its Moving Average
x = tsmovavg(x,'s',wind_size);
nan_ind = isnan(x);
x = x(~nan_ind);

%% ---- Fetting Matched Filter Output and Getting Its Moving Average  ----

y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;
y = y(g_ix); % Using matched filter outputs of good records.

% Calculating Its Moving Average
y = tsmovavg(y,'s',wind_size);
nan_ind = isnan(y);
y = y(~nan_ind);

%% ---- Performing Simple Linear Regression ----

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
xlabel('Moving Average of Temperature', 'FontSize', 14);
ylabel('Moving Average of Matched Filter Output', 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','NorthWest');
set(fig_leg,'FontSize',10);
hold off;
% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);
% ---- Saving Results ----
saveas(fig_id,[figDir,'/',card,'_',num2str(deg),'.pdf']);
% -------------------------------