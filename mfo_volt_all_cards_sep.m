% *** Mean (MFO) for All Cards ~ Mean (Supply Voltage) for All Cards ***
% This regression is performed not on the concatenation of datasets.
%% ---- Preparing Environment ----
clear all;
close all;

deg = 1; % Degree of Regression

tRec = 10000; % Number of Records (For The Same Datasets)

cards = char('b5c1','b5c2','b5c3','b5c5','b5c7');

% ---- Dataset 1 Variables ----
testDir1 = '/media/SHAYAN_HDD/Results/Collection_5/test_res_reg/1';
recsDir1 = '/media/SHAYAN_HDD/Data/collection_5';
sv_Tran1 = zeros(length(cards),tRec);
sv_IC1 = zeros(length(cards),tRec);
mfo_dat1 = zeros(length(cards),tRec);

Tran_mean1 = zeros(1,length(cards));
IC_mean1 = zeros(1,length(cards));
mfo_mean1 = zeros(1,length(cards));

% ---- Dataset 2 Variables ----
testDir2 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir2 = '/media/SHAYAN_HDD/Data/collection_6';
sv_Tran2 = zeros(length(cards),tRec);
sv_IC2 = zeros(length(cards),tRec);
mfo_dat2 = zeros(length(cards),tRec);

Tran_mean2 = zeros(1,length(cards));
IC_mean2 = zeros(1,length(cards));
mfo_mean2 = zeros(1,length(cards));

% ---- Directory of Results ----
figDir = '/media/SHAYAN_HDD/Results/Analysis_7/mfo_volt_all_cards_sep';

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data of Dataset 1 ----

disp('Getting Data for Dataset 1:');

for j = 1:length(cards)

    % ---- Getting Voltage Data ----
    for i = 1:tRec
    
        if mod(i,n) == 0
            disp([cards(j,:) ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
        end

        if (i >= 1) && (i <= 9)
            load([recsDir1,'/',cards(j,:),'/sample','0000',num2str(i),'.mat'],'volt_dat');
        elseif (i >= 10) && (i <= 99)
            load([recsDir1,'/',cards(j,:),'/sample','000',num2str(i),'.mat'],'volt_dat');
        elseif (i >= 100) && (i <= 999)
            load([recsDir1,'/',cards(j,:),'/sample','00',num2str(i),'.mat'],'volt_dat');
        elseif (i >= 1000) && (i <= 9999)
            load([recsDir1,'/',cards(j,:),'/sample','0',num2str(i),'.mat'],'volt_dat');
        else
            load([recsDir1,'/',cards(j,:),'/sample',num2str(i),'.mat'],'volt_dat');
        end
    
        if (isempty(volt_dat) == 0)
            sv_Tran1(j,i) = volt_dat(1);
            sv_IC1(j,i) = volt_dat(2);
        else
            sv_Tran1(j,i) = 0.0;
            sv_IC1(j,i) = 0.0;
        end
    
    end

    % Finding Good Records: The ones that have temperature data.
    [tran_temp,~] = outliers(sv_Tran1(j,:),10);
    nan_x = isnan(tran_temp);
    tran_temp = tran_temp(~nan_x);

    tran_ind1 = find (tran_temp~=0.0);
    Tran_mean1(j) = mean(tran_temp(tran_ind1));

    [ic_temp,~] = outliers(sv_IC1(j,:),10);
    nan_x = isnan(ic_temp);
    ic_temp = ic_temp(~nan_x);

    ic_ind1 = find (ic_temp~=0.0);
    IC_mean1(j) = mean(ic_temp(ic_ind1));

    % ---- Getting Matched Filter Output ----
    load([testDir1,'/',cards(j,:),'_',cards(j,:),'.mat'],'m0');
    mfo_dat1(j,:) = m0.genOrig;

    [mfo_temp,~] = outliers(mfo_dat1(j,:),10);
    nan_x = isnan(mfo_temp);
    mfo_temp = mfo_temp(~nan_x);

    mfo_mean1(j) = mean(mfo_temp);

end

%% ---- Getting Data of Dataset 2 ----

disp('Getting Data for Dataset 2:');

for j = 1:length(cards)

    % ---- Getting Voltage Data ----
    for i = 1:tRec
    
        if mod(i,n) == 0
            disp([cards(j,:) ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
        end

        if (i >= 1) && (i <= 9)
            load([recsDir2,'/',cards(j,:),'/sample','0000',num2str(i),'.mat'],'volt_dat');
        elseif (i >= 10) && (i <= 99)
            load([recsDir2,'/',cards(j,:),'/sample','000',num2str(i),'.mat'],'volt_dat');
        elseif (i >= 100) && (i <= 999)
            load([recsDir2,'/',cards(j,:),'/sample','00',num2str(i),'.mat'],'volt_dat');
        elseif (i >= 1000) && (i <= 9999)
            load([recsDir2,'/',cards(j,:),'/sample','0',num2str(i),'.mat'],'volt_dat');
        else
            load([recsDir2,'/',cards(j,:),'/sample',num2str(i),'.mat'],'volt_dat');
        end
    
        if (isempty(volt_dat) == 0)
            sv_Tran2(j,i) = volt_dat(1);
            sv_IC2(j,i) = volt_dat(2);
        else
            sv_Tran2(j,i) = 0.0;
            sv_IC2(j,i) = 0.0;
        end
    
    end

    % Finding Good Records: The ones that have temperature data.
    [tran_temp,~] = outliers(sv_Tran2(j,:),10);
    nan_x = isnan(tran_temp);
    tran_temp = tran_temp(~nan_x);

    tran_ind2 = find (tran_temp~=0.0);
    Tran_mean2(j) = mean(tran_temp(tran_ind2));

    [ic_temp,~] = outliers(sv_IC2(j,:),10);
    nan_x = isnan(ic_temp);
    ic_temp = ic_temp(~nan_x);

    ic_ind2 = find (ic_temp~=0.0);
    IC_mean2(j) = mean(ic_temp(ic_ind2));

    % ---- Getting Matched Filter Output ----
    load([testDir2,'/',cards(j,:),'_',cards(j,:),'.mat'],'m0');
    mfo_dat2(j,:) = m0.genOrig;

    [mfo_temp,~] = outliers(mfo_dat2(j,:),10);
    nan_x = isnan(mfo_temp);
    mfo_temp = mfo_temp(~nan_x);

    mfo_mean2(j) = mean(mfo_temp);

end

%% ---- Regression Analysis for Transformer Supply Voltage ----

% ---- Preparing Variables ----
x = [Tran_mean1,Tran_mean2];
y = [mfo_mean1,mfo_mean2];

% ---- Correlation Coefficient ----
cor_res = corrcoef(x',y');

% ---- Performing Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);

%% ---- Plotting Results ----

fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Mean of Transformer Supply Voltage (All Cards)', 'FontSize', 11);
ylabel('Mean of Matched Filter Output (All Cards)', 'FontSize', 11);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','NorthEast');
set(fig_leg,'FontSize',10);
hold off;
% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
title(['R^2 = ',Rsquare,' , ', 'Correlation Coefficient = ',num2str(cor_res(2))], 'FontSize', 12);
% ---- Saving Results ----
saveas(fig_id,[figDir,'/','all_cards','_','tran_',num2str(deg),'.pdf']);

%% ---- Regression Analysis for IC Supply Voltage ----

% ---- Preparing Variables ----
x = [IC_mean1,IC_mean2];
y = [mfo_mean1,mfo_mean2];

% ---- Correlation Coefficient ----
cor_res = corrcoef(x',y');

% ---- Performing Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);

%% ---- Plotting Results ----

fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Mean of IC Supply Voltage (All Cards)', 'FontSize', 11);
ylabel('Mean of Matched Filter Output (All Cards)', 'FontSize', 11);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','NorthEast');
set(fig_leg,'FontSize',10);
hold off;
% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
title(['R^2 = ',Rsquare,' , ', 'Correlation Coefficient = ',num2str(cor_res(2))], 'FontSize', 12);
% ---- Saving Results ----
saveas(fig_id,[figDir,'/','all_cards','_','ic_',num2str(deg),'.pdf']);
% -------------------------------