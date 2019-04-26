% *** Performing Multiple Linear Regression using Original Data ***
% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c7';
recsDir = '/media/SHAYAN_HDD/Collected_Signal/b5c7';
testDir = '/media/SHAYAN_HDD/Results/test_res/1';
tRec = 10000; % Number of Records
% ---- Fetching Temperature Data (First Predictor) ----
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

% ---- Fetching Actual Response Data  ----
y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
k = m0.genOrig; % It is used to obtain the second predictor.
y = m0.genOrig;
y = y(g_ix); % Using matched filter outputs of good records.

% --- Calculating Normalized Vector of Response Data (Second Predictor) ---
z = k/norm(k);
z = z(g_ix);
% ---- Preparing Data for Fitting ----
[x,z,y] = prepareSurfaceData(x,z,y);
% ---- Fitting Prepared Data ----
Fit_Object = fit([x,z],y,'poly11');
% ---- Getting Values of Coefficients ----
p_vals = coeffvalues(Fit_Object);
% ---- Calculating Estimated Response Data ----
yfit = p_vals(1) + p_vals(2)*x + p_vals(3)*z;
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
scatter3(x,z,y,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
xlabel('Temperature', 'FontSize', 14);
zlabel('Normalized Vector', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
set(gca,'YTickLabel',sprintf('%1.4f|',y));
title(['Multiple Linear Regression Analysis',' -- ','R^2 = ',num2str(rsq),' , ', 'R^2 Adjusted = ',num2str(rsq_adj)], 'FontSize', 14);
set(gca, 'fontsize', 12);
hold off;
% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/reg_res/',card,'.pdf']);
% -------------------------------