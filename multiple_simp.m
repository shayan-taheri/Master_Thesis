% *** Performing Multiple Linear Regression (MFO ~ Temp + Volt) ***
%% ---- Preparing Environment ----
clear all;
close all;

card = 'b5c7'; % Name of Card

tRec = 10000; % Number of Records

% ---- Dataset 1 Variables ----
testDir1 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir1 = ['/media/SHAYAN_HDD/Data/collection_6/',card];
temp_dat1 = zeros(1,tRec);
mfo_dat1 = zeros(1,tRec);
sv_Tran1 = zeros(1,tRec);
sv_IC1 = zeros(1,tRec);

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data of Dataset ----

disp('Getting Data for Dataset:');

% ---- Getting Temperature and Voltage Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'temp','volt_dat');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'temp','volt_dat');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'temp','volt_dat');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'temp','volt_dat');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'temp','volt_dat');
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
    
end

% Finding Good Records: The ones that have data.
temp_ind1 = find (temp_dat1~=0.0);
temp_dat1 = temp_dat1(temp_ind1);

sv_Tran1 = sv_Tran1(temp_ind1);

sv_IC1 = sv_IC1(temp_ind1);

% ---- Getting Matched Filter Output ----
load([testDir1,'/',card,'_',card,'.mat'],'m0');
mfo_dat1 = m0.genOrig;
mfo_dat1 = mfo_dat1(temp_ind1); % Using matched filter outputs of good records.

%% ---- Regression Analysis using Transformer Voltage ----

% ---- Preparing Data for Fitting ----
x = temp_dat1;
z1 = sv_Tran1;
y = mfo_dat1;
[x,z1,y] = prepareSurfaceData(x,z1,y);

% ---- Removing Outliers ----
[x,~] = outliers(x,10);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);
z1 = z1(~nan_x);

[y,~] = outliers(y,10);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
z1 = z1(~nan_y);

[z1,~] = outliers(z1,10);
nan_z1 = isnan(z1);
z1 = z1(~nan_z1);
y = y(~nan_z1);
x = x(~nan_z1);

% ---- Fitting Prepared Data ----
Fit_Object = fit([x,z1],y,'poly11');
% ---- Getting Values of Coefficients ----
p_vals = coeffvalues(Fit_Object);
% ---- Calculating Estimated Response Data ----
yfit = p_vals(1) + p_vals(2)*x + p_vals(3)*z1;
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p_vals));
% ---- Plotting Results ----
fig_id = figure();
plot(Fit_Object);
hold on;
scatter3(x,z1,y,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
xlabel('Temperature', 'FontSize', 14);
ylabel('Transformer Supply Voltage', 'FontSize', 14);
zlabel('Matched Filter Output', 'FontSize', 14);
title(['Multiple Linear Regression Analysis',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
hold off;
% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/multiple_simp/col6/',card,'_tran','.pdf']);
% ---- Plotting Residual ----
fig_id = figure();
stem(yresid);
xlabel('Number of Records', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title('Transformer Supply Voltage and Temperature', 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/multiple_simp/col6/',card,'_resi','_tran','.pdf']);

%% ---- Regression Analysis using IC Voltage ----

% ---- Preparing Data for Fitting ----
x = temp_dat1;
z2 = sv_IC1;
y = mfo_dat1;
[x,z2,y] = prepareSurfaceData(x,z2,y);

% ---- Removing Outliers ----
[x,~] = outliers(x,10);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);
z2 = z2(~nan_x);

[y,~] = outliers(y,10);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
z2 = z2(~nan_y);

[z2,~] = outliers(z2,10);
nan_z2 = isnan(z2);
z2 = z2(~nan_z2);
y = y(~nan_z2);
x = x(~nan_z2);

% ---- Fitting Prepared Data ----
Fit_Object = fit([x,z2],y,'poly11');
% ---- Getting Values of Coefficients ----
p_vals = coeffvalues(Fit_Object);
% ---- Calculating Estimated Response Data ----
yfit = p_vals(1) + p_vals(2)*x + p_vals(3)*z2;
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(p_vals));
% ---- Plotting Results ----
fig_id = figure();
plot(Fit_Object);
hold on;
scatter3(x,z2,y,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
xlabel('Temperature', 'FontSize', 14);
ylabel('IC Supply Voltage', 'FontSize', 14);
zlabel('Matched Filter Output', 'FontSize', 14);
title(['Multiple Linear Regression Analysis',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
hold off;
% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/multiple_simp/col6/',card,'_ic','.pdf']);
% ---- Plotting Residual ----
fig_id = figure();
stem(yresid);
xlabel('Number of Records', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title('IC Supply Voltage and Temperature', 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_3/multiple_simp/col6/',card,'_resi','_ic','.pdf']);
% -------------------------------