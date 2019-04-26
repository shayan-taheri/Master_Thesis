% *** Moving Average (MFO) ~ Moving Average (Supply Voltage) ***

% Correlation Coefficient (M.A.(MFO), M.A.(Supply Voltage))

%% ---- Preparing Environment ----
clear all;
close all;

card = 'b5c7';
deg = 1; % Degree of Regression
wind_size = 15; % Window Size for Moving Average

tRec = 10000; % Number of Records

% ---- Dataset 1 Variables ----
testDir1 = '/media/SHAYAN_HDD/Results/Collection_5/test_res_reg/1';
recsDir1 = ['/media/SHAYAN_HDD/Data/collection_5/',card];
mfo_dat1 = zeros(1,tRec);
sv_Tran1 = zeros(1,tRec);
sv_IC1 = zeros(1,tRec);

% ---- Dataset 2 Variables ----
testDir2 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir2 = ['/media/SHAYAN_HDD/Data/collection_6/',card];
mfo_dat2 = zeros(1,tRec);
sv_Tran2 = zeros(1,tRec);
sv_IC2 = zeros(1,tRec);

figDir = '/media/SHAYAN_HDD/Results/Analysis_6/MA_mfo_MA_volt';

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data of Dataset 1 ----

disp('Getting Data for Dataset 1:');

% ---- Getting Voltage Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'volt_dat');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'volt_dat');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'volt_dat');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'volt_dat');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'volt_dat');
    end
    
    if (isempty(volt_dat) == 0)
        sv_Tran1(i) = volt_dat(1);
        sv_IC1(i) = volt_dat(2);
    else
        sv_Tran1(i) = 0.0;
        sv_IC1(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
tran_ind1 = find (sv_Tran1~=0.0);
sv_Tran1 = sv_Tran1(tran_ind1);

ic_ind1 = find (sv_IC1~=0.0);
sv_IC1 = sv_IC1(ic_ind1);

% ---- Getting Matched Filter Output ----
load([testDir1,'/',card,'_',card,'.mat'],'m0');
mfo_dat1 = m0.genOrig;

%% ---- Getting Data of Dataset 2 ----

disp('Getting Data for Dataset 2:');

% ---- Getting Voltage Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir2,'/sample','0000',num2str(i),'.mat'],'volt_dat');
    elseif (i >= 10) && (i <= 99)
        load([recsDir2,'/sample','000',num2str(i),'.mat'],'volt_dat');
    elseif (i >= 100) && (i <= 999)
        load([recsDir2,'/sample','00',num2str(i),'.mat'],'volt_dat');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir2,'/sample','0',num2str(i),'.mat'],'volt_dat');
    else
        load([recsDir2,'/sample',num2str(i),'.mat'],'volt_dat');
    end
    
    if (isempty(volt_dat) == 0)
        sv_Tran2(i) = volt_dat(1);
        sv_IC2(i) = volt_dat(2);
    else
        sv_Tran2(i) = 0.0;
        sv_IC2(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
tran_ind2 = find (sv_Tran2~=0.0);
sv_Tran2 = sv_Tran2(tran_ind2);

ic_ind2 = find (sv_IC2~=0.0);
sv_IC2 = sv_IC2(ic_ind2);

% ---- Getting Matched Filter Output ----
load([testDir2,'/',card,'_',card,'.mat'],'m0');
mfo_dat2 = m0.genOrig;

%% ---- Regression Analysis for Transformer Supply Voltage ----

% ---- Preparing Variables ----
x = [sv_Tran1,sv_Tran2];
y = [mfo_dat1(tran_ind1),mfo_dat2(tran_ind2)];

% ---- Calculating Moving Average of Data ----
x = tsmovavg(x,'s',wind_size);
nan_ind = isnan(x);
x = x(~nan_ind);

y = tsmovavg(y,'s',wind_size);
nan_ind = isnan(y);
y = y(~nan_ind);

% ---- Removing Outliers ----
[x,~] = outliers(x,15);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,15);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

% ---- Correlation Coefficient ----
cor_res = corrcoef(x',y');

% ---- Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);

% ---- Plotting Results ----
fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Moving Average of Transformer Supply Voltage', 'FontSize', 12);
ylabel('Moving Average of Matched Filter Output', 'FontSize', 12);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
title(['Correlation Coefficient = ',num2str(cor_res(2)),' , ','R^2 = ',num2str(rsq),...
    ' (Degree ',num2str(deg),')'], 'FontSize', 12);

% ---- Saving Results ----
saveas(fig_id,[figDir,'/',card,'_','tran_',num2str(deg),'.pdf']);

%% ---- Regression Analysis for Supply Voltage of IC ----

% ---- Preparing Variables ----
x = [sv_IC1,sv_IC2];
y = [mfo_dat1(ic_ind1),mfo_dat2(ic_ind2)];

% ---- Calculating Moving Average of Data ----
x = tsmovavg(x,'s',wind_size);
nan_ind = isnan(x);
x = x(~nan_ind);

y = tsmovavg(y,'s',wind_size);
nan_ind = isnan(y);
y = y(~nan_ind);

% ---- Removing Outliers ----
[x,~] = outliers(x,15);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,15);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

% ---- Correlation Coefficient ----
cor_res = corrcoef(x',y');

% ---- Linear Regression ----
p = polyfit(x,y,deg);
yfit = polyval(p,x);

% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);

% ---- Plotting Results ----
fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Moving Average of IC Supply Voltage', 'FontSize', 12);
ylabel('Moving Average of Matched Filter Output', 'FontSize', 12);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
title(['Correlation Coefficient = ',num2str(cor_res(2)),' , ','R^2 = ',num2str(rsq),...
    ' (Degree ',num2str(deg),')'], 'FontSize', 12);

% ---- Saving Results ----
saveas(fig_id,[figDir,'/',card,'_','ic_',num2str(deg),'.pdf']);
% -------------------------------