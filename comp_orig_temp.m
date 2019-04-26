% *** Performing Comparison Between Datasets ***

%% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c7';

recsDir1 = ['/media/SHAYAN_HDD/Data/collection_1/',card];
testDir1 = '/media/SHAYAN_HDD/Results/Collection_1/test_res/1';

recsDir2 = ['/media/SHAYAN_HDD/Data/collection_2/',card];
testDir2 = '/media/SHAYAN_HDD/Results/Collection_2/test_res_reg/1';

tRec = 10000; % Number of Records
n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Fetching Data for Dataset 1 ----

disp('Fetching Data for Dataset 1:');

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

% ---- Getting Matched Filter Output ----
y1 = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir1,'/',card,'_',card,'.mat'],'m0');
y1 = m0.genOrig;
y1 = y1(g_ix1); % Using matched filter outputs of good records.

% ---- Correlation Coefficient ----
cor_res1 = corrcoef(x1',y1');

%% ---- Fetching Data for Dataset 2 ----

disp('Fetching Data for Dataset 2:');

% ---- Getting Temperature Data ----
x2 = zeros(1,tRec); % Vector for Temperature Data
for i = 1:tRec
    
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
    end

    if (i >= 1) && (i <= 9)
        load([recsDir2,'/sample','0000',num2str(i),'.mat'],'temp');
    elseif (i >= 10) && (i <= 99)
        load([recsDir2,'/sample','000',num2str(i),'.mat'],'temp');
    elseif (i >= 100) && (i <= 999)
        load([recsDir2,'/sample','00',num2str(i),'.mat'],'temp');
    elseif (i >= 1000) && (i <= 9999)
        load([recsDir2,'/sample','0',num2str(i),'.mat'],'temp');
    else
        load([recsDir2,'/sample',num2str(i),'.mat'],'temp');
    end
    
    if (isempty(temp) == 0)
        x2(i) = temp(1);
    else
        x2(i) = 0.0;
    end
    
end

% Finding Good Records: The ones that have temperature data.
g_ix2 = find (x2~=0.0);
x2 = x2(g_ix2);

% ---- Getting Matched Filter Output ----
y2 = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir2,'/',card,'_',card,'.mat'],'m0');
y2 = m0.genOrig;
y2 = y2(g_ix2); % Using matched filter outputs of good records.

% ---- Correlation Coefficient ----
cor_res2 = corrcoef(x2',y2');

%% ---- Plotting Overlapping Data ----
fig_id = figure();
plot(x1,y1,'*b');
hold on;
plot(x2,y2,'or');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title('Comparison - Overlapping','FontSize',18);
set(gca, 'fontsize', 12);
legend('Dataset 1','Dataset 2');
hold off;
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_2/reg_res/',card,'_','form1','.pdf']);

%% ---- Plotting Comparison and Correlation Coefficient ----
fig_id = figure();

subplot(2,1,1);
plot(x1,y1,'*b');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title('Dataset 1','FontSize',18);
set(gca, 'fontsize', 12);
text(0.30564,2.8278*10^5,['Correlation Coefficient:', ...
char(10),num2str(cor_res1(1)),'                 ',num2str(cor_res1(2)), ...
char(10),num2str(cor_res1(3)),'      ',num2str(cor_res1(4))], ...
'FontSize',12,'FontName','Times New Roman');

subplot(2,1,2);
plot(x2,y2,'or');
xlabel('Temperature', 'FontSize', 14);
ylabel('Matched Filter Output', 'FontSize', 14);
title('Dataset 2','FontSize',18);
set(gca, 'fontsize', 12);
text(0.31,2.824*10^5,['Correlation Coefficient:', ...
char(10),num2str(cor_res2(1)),'                 ',num2str(cor_res2(2)), ...
char(10),num2str(cor_res2(3)),'      ',num2str(cor_res2(4))], ...
'FontSize',12,'FontName','Times New Roman');
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_2/reg_res/',card,'_','form2','.pdf']);

% -------------------------------