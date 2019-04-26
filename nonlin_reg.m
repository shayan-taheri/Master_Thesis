% *** Performing Non-Linear Regression using Original Data ***
% Different Options for Model Function:
% 1) Exponential: b(1)*exp(b(2)*x(:,1))
% 2) Power: b(1)*(x(:,1).^b(2))
% 3) Saturation: (b(1)*x(:,1))./(b(2)+x(:,1))
% 4) Linear-Log: b(1)+b(2)*log(x(:,1))
% 5) Log-Linear: exp(b(1)+b(2)*x(:,1))
% Notice: They should be placed after @(b,x) statement.
% ---- Preparing Environment ----
clear all;
close all;
model_name = 'Saturation'; % Name of Model Function
modelfun = @(b,x) (b(1)*x(:,1))./(b(2)+x(:,1)); % Description of Function.
coef_num = 2; % Number of Coefficients.
card = 'b5c1';
recsDir = '/media/SHAYAN_HDD/Collected_Signal/b5c1';
testDir = '/media/SHAYAN_HDD/Results/test_res/1';
tRec = 10000; % Number of Records
beta0 = randn(coef_num,1); % Initial Values for Coefficient Vector
% ---- Fetching Temperature Data ----
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

% ---- Matched Filter Output ----
y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;
y = y(g_ix); % Using matched filter outputs of good records.
% ---- Preparing Data ----
x = x';
y = y';
% ---- Linear Regression ----
mdl = NonLinearModel.fit(x,y,modelfun,beta0);
yfit = feval(mdl,x);
% ---- R-Square Calculation ----
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - (SSresid/SStotal);
rsq_adj = 1 - (SSresid/SStotal)*(length(y)-1)/(length(y)-length(coef_num));
% ---- Plotting Results ----
fig_id = figure();
plot(x,y,'*b');
hold on;
plot(x,yfit,'-r');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title(['Non-Linear Regression Analysis - Model Function: ',model_name], 'FontSize', 14);
set(gca, 'fontsize', 12);
fig_leg = legend ('Observed Data','Estimated Data');
set(fig_leg,'FontSize',10);
hold off;
% ---- Representing R-Squared ----
Rsquare = num2str(rsq);
Rsquare_Adj = num2str(rsq_adj);
text(max(x)-0.011,min(y)+40,['R^2 = ',Rsquare,' , ', 'R^2 Adjusted = ',Rsquare_Adj],'FontSize',12,'FontName','Times New Roman');
% ---- Saving Results ----
saveas(fig_id,['/media/SHAYAN_HDD/Results/reg_res/',card,'_',model_name,'.pdf']);
% -------------------------------