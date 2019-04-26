% *** Thermal System Modeling (Bidirectional) + Simple Regression Analysis ***

% Equivalent IC Name: REALTEK RTL8019AS

% Actual IC Name: Genica GN-788

% ------------ Silicon ------------
% Thermal Conductivity = 149 Watt/(Meter * Kelvin)
% Specific Heat = 704.5985 Joule/(Kilogram * Kelvin)

% ------------ Plastic (Nylon 6,6) ------------
% Thermal Conductivity = 0.25 Watt/(Meter * Kelvin)
% Specific Heat = 1670 Joule/(Kilogram * Kelvin)

% ------------ IC Dimension and Weight ------------
% Length = 20 mm , Width = 14 mm , Height = 2.85 mm , Weight = 3.49 gram

% Thermal Conductivity Range = 0.25:0.9:148 --> 165 Numbers

% Specific Heat Range = 705:5.85:1670 --> 165 Numbers

%% ---- Preparing Environment ----
clear all;
close all;

deg = 3; % Degree of Regression

card = 'b5c7';

initial_temp = 0.2705556; % Initial Internal Temperature of Device

testDir = '/media/SHAYAN_HDD/Results/Collection_7/test_res_reg/1'; % Directory of MFO

recsDir1 = ['/media/SHAYAN_HDD_2/Data/2015-SPRING-07/',card]; % Directory of Records

resDir = '/media/SHAYAN_HDD/Results/TSM_bidirectional'; % Directory of Results

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

disp('Fetching Data for Dataset:');

% ---- Getting All Data ----

temp_vect = zeros(1,tRec); % Vector for Temperature Data
date_raw = zeros(1,tRec); % Vector for Date Data

for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'temp','ts');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'temp','ts');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'temp','ts');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'temp','ts');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'temp','ts');
    end
    
    date_raw(i) = ts;
    
    if (isempty(temp) == 0)
        temp_vect(i) = temp(1);
    else
        temp_vect(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
temp_ind = find (temp_vect~=0.0);
temp_vect = temp_vect(temp_ind);

% ---- Calculating Time Data ----
initial_time = 300; % Assumption for Start Time = 5 Minutes.

date_diff = diff(date_raw);

date_diff = datevec(date_diff);

date_diff(:,6) = date_diff(:,4) * 3600 + date_diff(:,5) * 60 + date_diff(:,6);

date_diff = date_diff(:,6);

date_diff = [initial_time,date_diff'];

time_dat = date_diff;

time_dat = time_dat(temp_ind);

% ---- Matched Filter Output ----
mfo = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo = m0.genOrig;
mfo = mfo(temp_ind);

%% ---- Thermal System Modeling ----

% ---- Required Parameters for Modeling ----

thermal_conduct = 0.25:0.9:148; % Thermal Conductivity
specific_heat = 705:5.85:1670; % Specific Heat

% Inverse relationship between Thermal Conductivity and Specific Heat.
specific_heat = fliplr(specific_heat);

length_dim = 20 * 10^(-3);
width = 14 * 10^(-3);
height = 2.85 * 10^(-3);
weight = 3.49 * 10^(-3);

internal_temp = zeros(size(specific_heat,2),size(temp_vect,2)+1);

internal_temp(:,1) = initial_temp;

% ---- Obtaining Internal Temperature of Device ----

for j = 1:size(specific_heat,2) % Number of elements in the vectors of thermal properties.
    
    Rt = height / (thermal_conduct(j) * width * length_dim);
    Ct = weight * specific_heat(j);
    
    for k = 1:size(temp_vect,2)
        
        % Internal Temperature Calculation
        internal_temp(j,k+1) = ((temp_vect(k)*(Rt*Ct+time_dat(k))) + (Rt*Ct*internal_temp(j,k)))/(time_dat(k) + 2*Rt*Ct);
    
    end
    
end

internal_temp = internal_temp(:,2:end); % Eliminating Initial Temperature from Matrix

%% ---- Regression Analysis on Both Ambient and Internal Temperature ----

% ---- Ambient Temperature ----
x = temp_vect;
y = mfo;

[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

p = polyfit(x,y,deg);
yfit_amb = polyval(p,x);

yresid = y - yfit_amb;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq_amb = 1 - (SSresid/SStotal);

% ---- Internal Temperature ----
yfit_int_all = zeros(size(internal_temp,1),size(mfo,2)-50);
rsq_int_all = zeros(1,size(internal_temp,1));

for m = 1:size(internal_temp,1)
    
    x = internal_temp(m,:);
    y = mfo;

    [x,~] = outliers(x,25);
    nan_x = isnan(x);
    x = x(~nan_x);
    y = y(~nan_x);

    [y,~] = outliers(y,25);
    nan_y = isnan(y);
    y = y(~nan_y);
    x = x(~nan_y);

    p = polyfit(x,y,deg);
    yfit_int = polyval(p,x);

    yresid = y - yfit_int;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq_int = 1 - (SSresid/SStotal);
    
    yfit_int_all(m,:) = yfit_int;
    rsq_int_all(m) = rsq_int;

end

tmpr_rsq = rsq_int_all; % Temporary R-Square Values
min_rsq = min(rsq_int_all); % Minimum R-Square Value

rsq_index = zeros(1,3); % Index of Best R-Square Values

for t = 1:3

    [~,index_temp] = max(tmpr_rsq);
    rsq_index(t) = index_temp;
    tmpr_rsq(index_temp) = min_rsq;

end

%% ---- Plotting Results of Regression Analyses ----

% ---- Ambient Temperature ----
x = temp_vect;
y = mfo;

[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);

[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);

fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit_amb,'-r');
xlabel('Ambient Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data','Location','NorthWest');
set(fig_leg,'FontSize',10);
hold off;
title(['Degree ',num2str(deg),': ','R^2 = ',num2str(rsq_amb)], 'FontSize', 14);
saveas(fig_id,[resDir,'/',card,'_','amb_','deg_',num2str(deg),'.pdf']);

% ---- Internal Temperature ----

for jj = 1:3

    x = internal_temp(rsq_index(jj),:);
    y = mfo;

    [x,~] = outliers(x,25);
    nan_x = isnan(x);
    x = x(~nan_x);
    y = y(~nan_x);

    [y,~] = outliers(y,25);
    nan_y = isnan(y);
    y = y(~nan_y);
    x = x(~nan_y);
    
    fig_id = figure();
    plot(x,y,'*b');
    hold on;
    plot(x,yfit_int_all(rsq_index(jj),:),'-r');
    xlabel('Internal Temperature', 'FontSize', 14);
    ylabel('Matched Filter Output', 'FontSize', 14);
    set(gca, 'fontsize', 12);
    fig_leg = legend ('Observed Data','Estimated Data','Location','NorthWest');
    set(fig_leg,'FontSize',10);
    hold off;
    title(['Degree ',num2str(deg),': ','R^2 = ',num2str(rsq_int_all(rsq_index(jj))),' , Choice = ',num2str(jj)], 'FontSize', 14);
    saveas(fig_id,[resDir,'/',card,'_','int_',num2str(jj),'_deg_',num2str(deg),'.pdf']);

end

%% ---- Storing Data ----

save([resDir,'/',card,'_deg_',num2str(deg),'_data.mat'], ...
        'initial_temp','temp_vect','mfo','internal_temp','thermal_conduct', ...
        'specific_heat','rsq_index','rsq_int_all','rsq_amb');

% -----------------------