% *** Performing Correlation Analysis (Temp ~ Volt) ***
%% ---- Preparing Environment ----
clear all;
close all;

card = 'b5c7'; % Name of Card

tRec = 10000; % Number of Records

% ---- Dataset Variables ----
recsDir = ['/media/SHAYAN_HDD/Data/collection_6/',card];
temp_dat = zeros(1,tRec);
sv_Tran = zeros(1,tRec);
sv_IC = zeros(1,tRec);

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data ----

disp('Getting Data:');

% ---- Getting Temperature and Voltage Data ----
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir,'/sample','0000',num2str(i),'.mat'],'temp','volt_dat');
    elseif (i >= 10) && (i <= 99)
        load([recsDir,'/sample','000',num2str(i),'.mat'],'temp','volt_dat');
    elseif (i >= 100) && (i <= 999)
        load([recsDir,'/sample','00',num2str(i),'.mat'],'temp','volt_dat');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir,'/sample','0',num2str(i),'.mat'],'temp','volt_dat');
    else
        load([recsDir,'/sample',num2str(i),'.mat'],'temp','volt_dat');
    end
    
    if (isempty(volt_dat) == 0)
        sv_Tran(i) = volt_dat(1);
        sv_IC(i) = volt_dat(2);
    else
        sv_Tran(i) = 0.0;
        sv_IC(i) = 0.0;
    end
    
    if (isempty(temp) == 0)
        temp_dat(i) = temp(1);
    else
        temp_dat(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have data.
temp_ind = find (temp_dat~=0.0);
temp_dat = temp_dat(temp_ind);

sv_Tran = sv_Tran(temp_ind);

sv_IC = sv_IC(temp_ind);

%% ---- Plotting Comparison and Correlation Coefficient ----

fig_id = figure();

subplot(2,1,1);
x = temp_dat;
y = sv_Tran;
[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);
[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
cor_res1 = corrcoef(x',y');
plot(x,y,'*b');
xlabel('Temperature', 'FontSize', 14);
ylabel('Transformer Supply Voltage', 'FontSize', 14);
title(['Comparison 1 - Correlation Coefficient = ',num2str(cor_res1(2))],'FontSize',14);
set(gca, 'fontsize', 12);

subplot(2,1,2);
x = temp_dat;
y = sv_IC;
[x,~] = outliers(x,25);
nan_x = isnan(x);
x = x(~nan_x);
y = y(~nan_x);
[y,~] = outliers(y,25);
nan_y = isnan(y);
y = y(~nan_y);
x = x(~nan_y);
cor_res2 = corrcoef(x',y');
plot(x,y,'or');
xlabel('Temperature', 'FontSize', 14);
ylabel('IC Supply Voltage', 'FontSize', 14);
title(['Comparison 2 - Correlation Coefficient = ',num2str(cor_res2(2))],'FontSize',14);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Analysis_2/corr_temp_volt/col6/',card,'.pdf']);
% ----------------------------------------------------------