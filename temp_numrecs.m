% *** Performing Comparison Between Temperature of Datasets ***

%% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c7';

recsDir1 = ['/media/SHAYAN_HDD/Data/collection_3/',card];

tRec = 5346; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset ----

disp('Fetching Data for Dataset:');

% ---- Getting Temperature Data ----
x1 = zeros(1,tRec); % Vector for Temperature Data
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir1,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir1,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir1,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir1,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir1,'/sample',num2str(i),'.mat'],'temp');
    end
    
    if (isempty(temp) == 0)
        x1(i) = temp(1);
    else
        x1(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
g_ix1 = find (x1~=0.0);
x1 = x1(g_ix1);

%% ---- Plotting Comparison ----
fig_id = figure();

plot(g_ix1,x1,'*b');
xlabel('Index of Record', 'FontSize', 14);
ylabel('Temperature', 'FontSize', 14);
title('Temperature Data','FontSize',18);
set(gca, 'fontsize', 12);

saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_3/temp/',card,'_','temp','.pdf']);

% -------------------------------