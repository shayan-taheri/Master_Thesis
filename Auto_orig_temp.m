% *** Performing Simple Linear Regression for N-Degree Polynomial using Original Data ***
% ---- Preparing Environment ----
clear all;
close all;
nom_temp = 25; % Nominal Temperature
order_var = 0.1; % Order of Variation for Nominal Temperature
card = 'b5c7';
deg = 1; % Degree of Regression
iter_num = 1000; % Number of Iteration for Analysis
testDir = '/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name';
recsDir = ['/media/SHAYAN_HDD/Data/collection_2/',card];
tRec = 10000; % Number of Records
% ---- Fetching Temperature Data ----
x = zeros(1,tRec); % Vector for Temperature Data
n = 1000; % For Updating the Status
nomStr = num2str(nom_temp);
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
act_temp = x(g_ix); % Vector of Actual Temperature

% ---- Matched Filter Output ----
y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;
y = y(g_ix); % Using matched filter outputs of good records.

rsq_vect = zeros(1,iter_num); % Vector of R-Squared Values

nom_vect = zeros(1,iter_num); % Vector of Nominal Temperature

for i1 = 1:iter_num
    
    % ---- Constructing Predictor ----
    x = (nom_temp-act_temp)/(nom_temp);
    
    % ---- Linear Regression ----
    p = polyfit(x,y,deg);
    yfit = polyval(p,x);
    
    % ---- R-Square Calculation ----
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - (SSresid/SStotal);
    
    % ---- Vector of R-Squared Values ----
    rsq_vect(i1) = rsq;
    
    % ---- Vector of Nominal Temperature ----
    nom_vect(i1) = nom_temp;
    
    % ---- Updating State of Nominal Temperature ----
    nom_temp = nom_temp + order_var;

end

% ---- Analysis of Results ----
[rsq_val,rsq_ind] = max(rsq_vect);
nom_max = nom_vect(rsq_ind);
mat_vect = [nom_vect;rsq_vect];

% ---- Plot for Maximum R-Squared Value ----
x = (nom_max-act_temp)/(nom_max);
p = polyfit(x,y,deg);
yfit = polyval(p,x);

fig_id = figure();
plot(x*(10^4),y,'*b');
hold on;
plot(x*(10^4),yfit,'-r');
xlabel('Predictor * 10^4', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title(['Simple Linear Regression Analysis - Degree ',num2str(deg)], 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data');
set(fig_leg,'FontSize',10);
hold off;
text(9968.5,2.825*10^5,['R^2 = ',num2str(rsq_val)],'FontSize',12,'FontName','Times New Roman');
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_2/auto_orig/',card,...
    '_',nomStr,'_',num2str(order_var),'.pdf']);

% ---- Saving All R-Squared Values and Nominal Temperatures ----
dash_line = '----------------------------------------------------------';

fileID = fopen(['/media/SHAYAN_HDD/Results/Collection_2/auto_orig/',card,...
    '_',nomStr,'_',num2str(order_var),'.txt'],'w');
fprintf(fileID,'%-25s \t %-s\r\n','Nominal Temperature','R-Squared Value');
fprintf(fileID,'%-s\r\n',dash_line);
fprintf(fileID,'%-25.15f \t %-25.15f\r\n',mat_vect);
fclose(fileID);

% -------------------------------