% *** Performing Simple Linear Regression (MFO ~ Temp) ***
%% ---- Preparing Environment ----
clear all;
close all;

deg = 1; % Degree of Regression
card = 'b5c7'; % Name of Card

tRec = 10000; % Number of Records (For The Same Datasets)

% ---- Dataset 1 Variables ----
testDir1 = '/media/SHAYAN_HDD/Results/Collection_1/test_res_mean/name';
recsDir1 = ['/home/shayan/collection_1/',card];
temp_dat1 = zeros(1,tRec);
mfo_dat1 = zeros(1,tRec);

% ---- Dataset 2 Variables ----
testDir2 = '/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name';
recsDir2 = ['/media/SHAYAN_HDD/Data/collection_2/',card];
temp_dat2 = zeros(1,tRec);
mfo_dat2 = zeros(1,tRec);

% ---- Dataset 5 Variables ----
testDir3 = '/media/SHAYAN_HDD/Results/Collection_5/test_res_reg/1';
recsDir3 = ['/media/SHAYAN_HDD/Data/collection_5/',card];
temp_dat3 = zeros(1,tRec);
mfo_dat3 = zeros(1,tRec);

% ---- Dataset 6 Variables ----
testDir4 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir4 = ['/media/SHAYAN_HDD/Data/collection_6/',card];
temp_dat4 = zeros(1,tRec);
mfo_dat4 = zeros(1,tRec);

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data of Dataset 1 ----

disp('Getting Data for Dataset 1:');

% ---- Getting Temperature Data ----
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
        temp_dat1(i) = temp(1);
    else
        temp_dat1(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
temp_ind1 = find (temp_dat1~=0.0);
temp_dat1 = temp_dat1(temp_ind1);

% ---- Getting Matched Filter Output ----
load([testDir1,'/',card,'_',card,'.mat'],'m0');
mfo_dat1 = m0.genOrig;
mfo_dat1 = mfo_dat1(temp_ind1); % Using matched filter outputs of good records.

%% ---- Getting Data of Dataset 2 ----

disp('Getting Data for Dataset 2:');

% ---- Getting Temperature Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir2,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir2,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir2,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir2,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir2,'/sample',num2str(i),'.mat'],'temp');
    end
    
    if (isempty(temp) == 0)
        temp_dat2(i) = temp(1);
    else
        temp_dat2(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
temp_ind2 = find (temp_dat2~=0.0);
temp_dat2 = temp_dat2(temp_ind2);

% ---- Getting Matched Filter Output ----
load([testDir2,'/',card,'_',card,'.mat'],'m0');
mfo_dat2 = m0.genOrig;
mfo_dat2 = mfo_dat2(temp_ind2); % Using matched filter outputs of good records.

%% ---- Getting Data of Dataset 3 ----

disp('Getting Data for Dataset 3:');

% ---- Getting Temperature Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir3,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir3,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir3,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir3,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir3,'/sample',num2str(i),'.mat'],'temp');
    end
    
    if (isempty(temp) == 0)
        temp_dat3(i) = temp(1);
    else
        temp_dat3(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
temp_ind3 = find (temp_dat3~=0.0);
temp_dat3 = temp_dat3(temp_ind3);

% ---- Getting Matched Filter Output ----
load([testDir3,'/',card,'_',card,'.mat'],'m0');
mfo_dat3 = m0.genOrig;
mfo_dat3 = mfo_dat3(temp_ind3); % Using matched filter outputs of good records.

%% ---- Getting Data of Dataset 4 ----

disp('Getting Data for Dataset 4:');

% ---- Getting Temperature Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir4,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir4,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir4,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir4,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir4,'/sample',num2str(i),'.mat'],'temp');
    end
    
    if (isempty(temp) == 0)
        temp_dat4(i) = temp(1);
    else
        temp_dat4(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
temp_ind4 = find (temp_dat4~=0.0);
temp_dat4 = temp_dat4(temp_ind4);

% ---- Getting Matched Filter Output ----
load([testDir4,'/',card,'_',card,'.mat'],'m0');
mfo_dat4 = m0.genOrig;
mfo_dat4 = mfo_dat4(temp_ind4); % Using matched filter outputs of good records.

%% ---- Regression Analysis ----

% ---- Preparing Variables ----
x = [temp_dat1,temp_dat2,temp_dat3,temp_dat4];
y = [mfo_dat1,mfo_dat2,mfo_dat3,mfo_dat4];

% ---- Removing Outliers ----
[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

% ---- Performing Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p));

%% ---- Plotting Results ----

fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Temperature for Four Datasets', 'FontSize', 11);
ylabel('Matched Filter Output for Four Datasets', 'FontSize', 11);
title(['Simple Linear Regression Analysis - Degree ',num2str(deg)], 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthWest');
set(fig_leg,'FontSize',10);
hold off;
% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);
% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_2/Orig4_Temp4/',card,'_',num2str(deg),'.pdf']);
% -------------------------------