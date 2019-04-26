% *** Plotting Supply Voltages of Datasets ***

%% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c7';

recsDir1 = ['/media/SHAYAN_HDD/Data/collection_5/',card];

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

disp('Fetching Data for Dataset:');

% ---- Getting Supply Voltage Data ----

sv_Tran = zeros(1,tRec); % Vector for Supply Voltage of Transformer
sv_IC = zeros(1,tRec); % Vector for Supply Voltage of IC

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
        sv_Tran(i) = volt_dat(1);
        sv_IC(i) = volt_dat(2);
    else
        sv_Tran(i) = 0.0;
        sv_IC(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have supply voltge data.
tran_ind = find (sv_Tran~=0.0);
sv_Tran = sv_Tran(tran_ind);

ic_ind = find (sv_IC~=0.0);
sv_IC = sv_IC(ic_ind);

%% ---- Plotting Comparison ----
fig_id = figure();

subplot(2,1,1);
plot(tran_ind,sv_Tran,'*b');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Supply Voltage (V)', 'FontSize', 14);
title('Transformer Supply Voltage Analysis','FontSize',18);
set(gca, 'fontsize', 12);

subplot(2,1,2);
plot(ic_ind,sv_IC,'*r');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Supply Voltage (V)', 'FontSize', 14);
title('IC Supply Voltage Analysis','FontSize',18);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_5/voltage/',card,'.pdf']);

% -------------------------------