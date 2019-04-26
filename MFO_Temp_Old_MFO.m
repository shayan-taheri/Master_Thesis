% *** MFO (Now) ~ Temperature (Now) + MFO (Old) ***

%% ---- Preparing Environment ----

clear all;
close all;

card = 'b5c7';

% ---- Current MFO and Temperature Data ----
testDir = '/media/SHAYAN_HDD/Results/Collection_8/test_res_reg/1';
recsDir = ['/media/SHAYAN_HDD/Data/2015-SPRING-08/',card];

% ---- Previously Taken MFO Data ----
testDir_old = '/media/SHAYAN_HDD/Results/Collection_7/test_res_reg/1';

% ---- Result Directory ----
resDir = '/media/SHAYAN_HDD/Results/Analysis_14/MFO_Temp_Old_MFO/col8';

tRec = 10000; % Number of Records

%% ---- Fetching Current Temperature Data ----

temp_vect = zeros(1,tRec); % Vector for Current Temperature Data

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
        temp_vect(i) = temp(1);
    else
        temp_vect(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
temp_ind = find(temp_vect~=0.0);
temp_vect = temp_vect(temp_ind);
x = temp_vect;

%% ---- Fetching Current Matched Filter Output Data  ----

y = zeros(1,tRec); % Vector for Current Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;
y = y(temp_ind); % Using matched filter outputs of good records.

%% ---- Fetching Old Matched Filter Output Data ----

z = zeros(1,tRec); % Vector for Old Matched Filter Output
load([testDir_old,'/',card,'_',card,'.mat'],'m0');
z = m0.genOrig;
z = z(temp_ind); % Using matched filter outputs of good records.

%% ---- Performing Multiple Linear Regression ----

% ---- Removing Outliers ----
[x,~] = outliers(x,30); % Outliers of Temperature
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);
z = z(~nan_x);

[y,~] = outliers(y,30); % Outliers of Current MFO
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
z = z(~nan_y);

[z,~] = outliers(z,30); % Outliers of Old MFO
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
xlabel('Temperature (Now)', 'FontSize', 11);
ylabel('Matched Filter Output (Old)', 'FontSize', 11);
zlabel('Matched Filter Output (Now)', 'FontSize', 11);
title(['Multiple Linear Regression Analysis',' -- ','R^2 = ',num2str(rsq)], 'FontSize', 14);
set(gca, 'fontsize', 12);
hold off;
% ---- Saving Results ----
saveas(fig_id,[resDir,'/',card,'.pdf']);
% ---- Plotting Residual ----
fig_id = figure();
plot(yresid,'*');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Residual of Response Data', 'FontSize', 14);
title('Analysis of Residuals', 'FontSize', 14);
set(gca, 'fontsize', 12);
% ---- Saving Residual ----
saveas(fig_id,[resDir,'/',card,'_res','.pdf']);
% -------------------------------