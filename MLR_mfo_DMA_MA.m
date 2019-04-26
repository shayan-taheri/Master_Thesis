% *** MFO ~ Deriv. (M.A. (Temperature)) + Moving Average (Temperature) ***

%% ---- Preparing Environment ----

clear all;
close all;

deg = 1; % Degree of Regression
wind_size = 15; % Window Size for Moving Average

card = 'b5c7';

testDir = '/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name';
recsDir = ['/media/SHAYAN_HDD/Data/collection_2/',card];
figDir = '/media/SHAYAN_HDD/Results/Analysis_6/MLR_mfo_DMA_MA/col2';

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

z = diff(x); % Derivation of the Moving Average of Temperature.

x = x(2:end); % Trimming temperature data for derivation.

%% ---- Fetting Matched Filter Output and Getting Its Moving Average  ----

y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;
y = y(g_ix); % Using matched filter outputs of good records.

% Trimming data to be matched to the regressor data.
y = y(~nan_ind); % For moving average.
y = y(2:end); % For derivation.

%% ---- Performing Multiple Linear Regression ----

% ---- Removing Outliers ----
[x,~] = outliers(x,10);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);
z = z(~nan_x);

[y,~] = outliers(y,10);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
z = z(~nan_y);

[z,~] = outliers(z,10);
nan_z = isnan(z);
z = z(~nan_z);
y = y(~nan_z);
x = x(~nan_z);

% ---- Linear Regression ----
[x,z,y] = prepareSurfaceData(x,z,y);
Fit_Object = fit([x,z],y,'poly11');
p_vals = coeffvalues(Fit_Object);
yfit = p_vals(1) + p_vals(2)*x + p_vals(3)*z;
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
% ---- Plotting Results ----
fig_id = figure();
plot(Fit_Object);
hold on;
scatter3(x,z,y,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
xlabel('M.A. of Temperature', 'FontSize', 11);
ylabel('Deriv. of M.A. of Temperature', 'FontSize', 11);
zlabel('Matched Filter Output', 'FontSize', 11);
title(['Multiple Linear Regression Analysis',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
hold off;
% ---- Saving Results ----
saveas(fig_id,[figDir,'/',card,'.pdf']);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title('Analysis of Residuals', 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[figDir,'/',card,'_resi','.pdf']);
% -------------------------------