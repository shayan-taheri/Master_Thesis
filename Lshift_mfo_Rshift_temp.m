% *** Left Shifted MFO ~ Right Shifted Temperature ***

%% ---- Preparing Environment ----
clear all;
close all;

num_shift = 15; % Number of Elements to Shift
deg = 1; % Degree of Regression
card = 'b5c7';

testDir = '/media/SHAYAN_HDD/Results/Collection_1/test_res_mean/name';
recsDir1 = ['/home/shayan/collection_1/',card];
figDir = '/media/SHAYAN_HDD/Results/Analysis_6/Lshift_mfo_Rshift_temp/col1';

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

disp('Fetching Data for Dataset:');

% ---- Getting All Data ----

temp_vect = zeros(1,tRec); % Vector for Temperature Data

for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'temp');
    end

    if (isempty(temp) == 0)
        temp_vect(i) = temp(1);
    else
        temp_vect(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
temp_ind = find (temp_vect~=0.0);
temp_vect = temp_vect(temp_ind);

% ---- Matched Filter Output ----
mfo = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo = m0.genOrig;

%% ---- Performing Regression Analysis ----

% ---- Preparing Variables ----
x = temp_vect;
y = mfo(temp_ind); % Using matched filter outputs of good records.

% ---- Shifting Data ----
x = x(1:end-num_shift);
y = y(num_shift+1:end);

% ---- Removing Outliers ----
[x,~] = outliers(x,10);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,10);
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
xlabel('Shifted Temperature', 'FontSize', 14);
ylabel('Shifted Matched Filter Output', 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,[figDir,'/',card,'_',num2str(deg),'.pdf']);
% -------------------------------