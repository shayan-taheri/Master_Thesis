% *** Storing All Data (MFO, Temperature, and Voltage) + Mean (Data) ***
%% ---- Preparing Environment ----
clear all;
close all;

tRec = 10000; % Number of Records (For The Same Datasets)

cards = char('b5c1','b5c2','b5c3','b5c5','b5c7');

% ---- Dataset 1 Variables ----
testDir1 = '/media/SHAYAN_HDD/Results/Collection_5/test_res_reg/1';
recsDir1 = '/media/SHAYAN_HDD/Data/collection_5';
sv_Tran1 = zeros(1,tRec);
sv_IC1 = zeros(1,tRec);
mfo_dat1 = zeros(1,tRec);
temp_dat1 = zeros(1,tRec);

Tran_mean1 = zeros(1,length(cards));
IC_mean1 = zeros(1,length(cards));
mfo_mean1 = zeros(1,length(cards));
temp_mean1 = zeros(1,length(cards));

% ---- Dataset 2 Variables ----
testDir2 = '/media/SHAYAN_HDD/Results/Collection_6/test_res_reg';
recsDir2 = '/media/SHAYAN_HDD/Data/collection_6';
sv_Tran2 = zeros(1,tRec);
sv_IC2 = zeros(1,tRec);
mfo_dat2 = zeros(1,tRec);
temp_dat2 = zeros(1,tRec);

Tran_mean2 = zeros(1,length(cards));
IC_mean2 = zeros(1,length(cards));
mfo_mean2 = zeros(1,length(cards));
temp_mean2 = zeros(1,length(cards));

% ---- Directory of Results ----
ResDir = '/media/SHAYAN_HDD/Results/Analysis_8';

n = 1000; % For Updating the Status
tRecStr = num2str(tRec);

%% ---- Getting Data of Dataset 1 ----

disp('Getting Data for Dataset 1:');

for j = 1:length(cards)

    % ---- Getting Temperature and Voltage Data ----
    for i = 1:tRec
    
        if mod(i,n) == 0
            disp([cards(j,:) ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
        end

        if (i >= 1) && (i <= 9)
            load([recsDir1,'/',cards(j,:),'/sample','0000',num2str(i),'.mat'],'volt_dat','temp');
        elseif (i >= 10) && (i <= 99)
            load([recsDir1,'/',cards(j,:),'/sample','000',num2str(i),'.mat'],'volt_dat','temp');
        elseif (i >= 100) && (i <= 999)
            load([recsDir1,'/',cards(j,:),'/sample','00',num2str(i),'.mat'],'volt_dat','temp');
        elseif (i >= 1000) && (i <= 9999)
            load([recsDir1,'/',cards(j,:),'/sample','0',num2str(i),'.mat'],'volt_dat','temp');
        else
            load([recsDir1,'/',cards(j,:),'/sample',num2str(i),'.mat'],'volt_dat','temp');
        end
        
        if (isempty(volt_dat) == 0)
            sv_Tran1(i) = volt_dat(1);
            sv_IC1(i) = volt_dat(2);
        else
            sv_Tran1(i) = 0.0;
            sv_IC1(i) = 0.0;
        end
    
        if (isempty(temp) == 0)
            temp_dat1(i) = temp(1);
        else
            temp_dat1(i) = 0.0;
        end
        
    end
    
    % ---- Getting Matched Filter Output ----
    load([testDir1,'/',cards(j,:),'_',cards(j,:),'.mat'],'m0');
    mfo_dat1 = m0.genOrig;

    % ---- Finding Available Records of Temperature Data ----
    temp_ind = find(temp_dat1~=0.0);
    
    temp_dat1 = temp_dat1(temp_ind);
    sv_Tran1 = sv_Tran1(temp_ind);
    sv_IC1 = sv_IC1(temp_ind);
    mfo_dat1 = mfo_dat1(temp_ind);
    
    % ---- Removing Outlier (MFO) ----
    [mfo_dat1,~] = outliers(mfo_dat1,5);
    nan_x = isnan(mfo_dat1);
    
    temp_dat1 = temp_dat1(~nan_x);
    sv_Tran1 = sv_Tran1(~nan_x);
    sv_IC1 = sv_IC1(~nan_x);
    mfo_dat1 = mfo_dat1(~nan_x);
    
    % ---- Removing Outlier (Temperature) ----
    [temp_dat1,~] = outliers(temp_dat1,5);
    nan_x = isnan(temp_dat1);
    
    temp_dat1 = temp_dat1(~nan_x);
    sv_Tran1 = sv_Tran1(~nan_x);
    sv_IC1 = sv_IC1(~nan_x);
    mfo_dat1 = mfo_dat1(~nan_x);
    
    % ---- Removing Outlier (IC Supply Voltage) ----
    [sv_IC1,~] = outliers(sv_IC1,5);
    nan_x = isnan(sv_IC1);
    
    temp_dat1 = temp_dat1(~nan_x);
    sv_Tran1 = sv_Tran1(~nan_x);
    sv_IC1 = sv_IC1(~nan_x);
    mfo_dat1 = mfo_dat1(~nan_x);

    % ---- Removing Outlier (Transformer Supply Voltage) ----
    [sv_Tran1,~] = outliers(sv_Tran1,5);
    nan_x = isnan(sv_Tran1);
    
    temp_dat1 = temp_dat1(~nan_x);
    sv_Tran1 = sv_Tran1(~nan_x);
    sv_IC1 = sv_IC1(~nan_x);
    mfo_dat1 = mfo_dat1(~nan_x);
    
    % ---- Obtaining Mean of Data ----
    Tran_mean1(j) = mean(sv_Tran1);
    IC_mean1(j) = mean(sv_IC1);
    mfo_mean1(j) = mean(mfo_dat1);
    temp_mean1(j) = mean(temp_dat1);
    
	% ---- Storing Data ----
    Temp_data = temp_dat1;
    Tran_volt = sv_Tran1;
    IC_volt = sv_IC1;
    MFO_data = mfo_dat1;
    
    save([ResDir,'/','set_5_',cards(j,:),'_data.mat'], ...
        'Temp_data','Tran_volt','IC_volt','MFO_data');

end

%% ---- Getting Data of Dataset 2 ----

disp('Getting Data for Dataset 2:');

for j = 1:length(cards)

    % ---- Getting Temperature and Voltage Data ----
    for i = 1:tRec
    
        if mod(i,n) == 0
            disp([cards(j,:) ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs analyzed)...']);
        end

        if (i >= 1) && (i <= 9)
            load([recsDir2,'/',cards(j,:),'/sample','0000',num2str(i),'.mat'],'volt_dat','temp');
        elseif (i >= 10) && (i <= 99)
            load([recsDir2,'/',cards(j,:),'/sample','000',num2str(i),'.mat'],'volt_dat','temp');
        elseif (i >= 100) && (i <= 999)
            load([recsDir2,'/',cards(j,:),'/sample','00',num2str(i),'.mat'],'volt_dat','temp');
        elseif (i >= 1000) && (i <= 9999)
            load([recsDir2,'/',cards(j,:),'/sample','0',num2str(i),'.mat'],'volt_dat','temp');
        else
            load([recsDir2,'/',cards(j,:),'/sample',num2str(i),'.mat'],'volt_dat','temp');
        end
        
        if (isempty(volt_dat) == 0)
            sv_Tran2(i) = volt_dat(1);
            sv_IC2(i) = volt_dat(2);
        else
            sv_Tran2(i) = 0.0;
            sv_IC2(i) = 0.0;
        end
    
        if (isempty(temp) == 0)
            temp_dat2(i) = temp(1);
        else
            temp_dat2(i) = 0.0;
        end
        
    end
    
    % ---- Getting Matched Filter Output ----
    load([testDir2,'/',cards(j,:),'_',cards(j,:),'.mat'],'m0');
    mfo_dat2 = m0.genOrig;

    % ---- Finding Available Records of Temperature Data ----
    temp_ind = find(temp_dat2~=0.0);
    
    temp_dat2 = temp_dat2(temp_ind);
    sv_Tran2 = sv_Tran2(temp_ind);
    sv_IC2 = sv_IC2(temp_ind);
    mfo_dat2 = mfo_dat2(temp_ind);
    
    % ---- Removing Outlier (MFO) ----
    [mfo_dat2,~] = outliers(mfo_dat2,5);
    nan_x = isnan(mfo_dat2);
    
    temp_dat2 = temp_dat2(~nan_x);
    sv_Tran2 = sv_Tran2(~nan_x);
    sv_IC2 = sv_IC2(~nan_x);
    mfo_dat2 = mfo_dat2(~nan_x);
    
    % ---- Removing Outlier (Temperature) ----
    [temp_dat2,~] = outliers(temp_dat2,5);
    nan_x = isnan(temp_dat2);
    
    temp_dat2 = temp_dat2(~nan_x);
    sv_Tran2 = sv_Tran2(~nan_x);
    sv_IC2 = sv_IC2(~nan_x);
    mfo_dat2 = mfo_dat2(~nan_x);
    
    % ---- Removing Outlier (IC Supply Voltage) ----
    [sv_IC2,~] = outliers(sv_IC2,5);
    nan_x = isnan(sv_IC2);
    
    temp_dat2 = temp_dat2(~nan_x);
    sv_Tran2 = sv_Tran2(~nan_x);
    sv_IC2 = sv_IC2(~nan_x);
    mfo_dat2 = mfo_dat2(~nan_x);

    % ---- Removing Outlier (Transformer Supply Voltage) ----
    [sv_Tran2,~] = outliers(sv_Tran2,5);
    nan_x = isnan(sv_Tran2);
    
    temp_dat2 = temp_dat2(~nan_x);
    sv_Tran2 = sv_Tran2(~nan_x);
    sv_IC2 = sv_IC2(~nan_x);
    mfo_dat2 = mfo_dat2(~nan_x);
    
    % ---- Obtaining Mean of Data ----
    Tran_mean2(j) = mean(sv_Tran2);
    IC_mean2(j) = mean(sv_IC2);
    mfo_mean2(j) = mean(mfo_dat2);
    temp_mean2(j) = mean(temp_dat2);
    
	% ---- Storing Data ----
    Temp_data = temp_dat2;
    Tran_volt = sv_Tran2;
    IC_volt = sv_IC2;
    MFO_data = mfo_dat2;
    save([ResDir,'/','set_6_',cards(j,:),'_data.mat'], ...
        'Temp_data','Tran_volt','IC_volt','MFO_data');

end

%% ---- Storing All Mean Data ----
Tran_mean = [Tran_mean1',Tran_mean2'];
IC_mean = [IC_mean1',IC_mean2'];
MFO_mean = [mfo_mean1',mfo_mean2'];
Temp_mean = [temp_mean1',temp_mean2'];

save([ResDir,'/','mean_all_data.mat'], ...
    'Tran_mean','IC_mean','MFO_mean','Temp_mean');
% -------------------------------