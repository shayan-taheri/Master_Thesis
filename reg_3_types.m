% *** Performing Two Regressions (Temperature , Supply Voltages) ***
% Type 1) MFO ~ (Temperature / Supply Voltage)
% Type 2) MFO ~ (Temperature * Supply Voltage)

%% ---- Preparing Environment ----
clear all;
close all;
deg = 1; % Degree of Regression
card = 'b5c7';

testDir = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';

recsDir1 = ['/media/SHAYAN_HDD/Data/collection_6/',card];

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

disp('Fetching Data for Dataset:');

% ---- Getting All Data ----

sv_Tran = zeros(1,tRec); % Vector for Supply Voltage of Transformer
sv_IC = zeros(1,tRec); % Vector for Supply Voltage of IC
temp_vect = zeros(1,tRec); % Vector for Temperature Data

for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'volt_dat','temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'volt_dat','temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'volt_dat','temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'volt_dat','temp');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'volt_dat','temp');
    end
    
    if (isempty(volt_dat) == 0)
        sv_Tran(i) = volt_dat(1);
        sv_IC(i) = volt_dat(2);
    else
        sv_Tran(i) = 0.0;
        sv_IC(i) = 0.0;
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

sv_Tran = sv_Tran(temp_ind);

sv_IC = sv_IC(temp_ind);

% ---- Matched Filter Output ----
mfo = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo = m0.genOrig;

mfo = mfo(temp_ind);

%% ---- Regressions for Transformer Supply Voltage ----

% ---- Regression Type 1 ----

% ---- Preparing Variables ----
x = temp_vect ./ sv_Tran;
y = mfo;

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
xlabel('Temperature / Transformer Voltage', 'FontSize', 12);
ylabel('Matched Filter Output', 'FontSize', 12);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/col6/reg_tran/',card,'_','type1_',num2str(deg),'.pdf']);

% ---- Regression Type 2 ----

% ---- Preparing Variables ----
x = temp_vect .* sv_Tran;
y = mfo;

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
xlabel('Temperature * Transformer Voltage', 'FontSize', 12);
ylabel('Matched Filter Output', 'FontSize', 12);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/col6/reg_tran/',card,'_','type2_',num2str(deg),'.pdf']);

%% ---- Regressions for IC Supply Voltage ----

% ---- Regression Type 1 ----

% ---- Preparing Variables ----
x = temp_vect ./ sv_IC;
y = mfo;

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
xlabel('Temperature / IC Voltage', 'FontSize', 12);
ylabel('Matched Filter Output', 'FontSize', 12);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/col6/reg_ic/',card,'_','type1_',num2str(deg),'.pdf']);

% ---- Regression Type 2 ----

% ---- Preparing Variables ----
x = temp_vect .* sv_IC;
y = mfo;

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
xlabel('Temperature * IC Voltage', 'FontSize', 12);
ylabel('Matched Filter Output', 'FontSize', 12);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','SouthEast');
set(fig_leg,'FontSize',10);
hold off;

% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
title(['Degree ',num2str(deg),': ','R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj], 'FontSize', 14);

% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/col6/reg_ic/',card,'_','type2_',num2str(deg),'.pdf']);
% ------------------------