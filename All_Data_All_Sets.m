% *** Storing All Data (MFO, Temperature, Noisy Part) + Voltage ***
%% ---- Preparing Environment ----
clear all;
close all;

tRec1 = 10000; % Number of Records (Datasets Type 1 and 3)

tRec2 = 5305; % Number of Records (Datasets Type 2)

tSamp = 37900; % Number of Noisy Sample Points

cards = char('b5c1','b5c2','b5c3','b5c5','b5c7');

% ************** IMPORTANT **************
% The Length Of Directory Members Of A Directory Variable ...
% Should Be The Same.

% ---- Dataset Type 1: 10,000 Records without Supply Voltage ----
testDir1 = char('/media/SHAYAN_HDD/Results/Collection_1/test_res_mean/name','/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name');
recsDir1 = char('/media/SHAYAN_HDD_2/Data/2015-SPRING-01','/media/SHAYAN_HDD_2/Data/2015-SPRING-02');

% ---- Dataset Type 2: ~5,000 Records without Supply Voltage ----
testDir2 = char('/media/SHAYAN_HDD/Results/Collection_3/test_res_reg/1','/media/SHAYAN_HDD/Results/Collection_4/test_res_reg/1');
recsDir2 = char('/media/SHAYAN_HDD_2/Data/2015-SPRING-03','/media/SHAYAN_HDD_2/Data/2015-SPRING-04');

% ---- Dataset Type 3: 10,000 Records with Supply Voltage ----
testDir3 = char('/media/SHAYAN_HDD/Results/Collection_5/test_res_reg/1','/media/SHAYAN_HDD/Results/Collection_6/test_res_reg/1');
recsDir3 = char('/media/SHAYAN_HDD_2/Data/2015-SPRING-05','/media/SHAYAN_HDD_2/Data/2015-SPRING-06');

% ---- Required Variables for Data Extraction (Type 1 and 3) ----
mfo_dat1 = zeros(1,tRec1);
temp_dat1 = zeros(1,tRec1);
noisy_mean1 = zeros(1,tRec1);
noisy_std1 = zeros(1,tRec1);
sv_Tran1 = zeros(1,tRec1);
sv_IC1 = zeros(1,tRec1);

n1 = 1000; % For Updating the Status
tRecStr1 = num2str(tRec1);

% ---- Required Variables for Data Extraction (Type 2) ----
mfo_dat2 = zeros(1,tRec2);
temp_dat2 = zeros(1,tRec2);
noisy_mean2 = zeros(1,tRec2);
noisy_std2 = zeros(1,tRec2);

n2 = 1000; % For Updating the Status
tRecStr2 = num2str(tRec2);

% ---- Directory of Results ----
ResDir = '/media/SHAYAN_HDD/Results/Analysis_9';

%% ---- Getting Data of Dataset (Type 1) ----

disp('Getting Data for Dataset (Type 1):');

for k = 1:size(recsDir1,1)

    for j = 1:size(cards,1)

        % ---- Getting Temperature Data ----
        for i = 1:tRec1
    
            if mod(i,n1) == 0
                disp([cards(j,:) ' ' num2str(i/tRec1*100) '% complete (' num2str(i) '/' tRecStr1 ' recs analyzed)...']);
            end

            if (i >= 1) && (i <= 9)
                load([recsDir1(k,:),'/',cards(j,:),'/sample','0000',num2str(i),'.mat'],'temp','rec');
            elseif (i >= 10) && (i <= 99)
                load([recsDir1(k,:),'/',cards(j,:),'/sample','000',num2str(i),'.mat'],'temp','rec');
            elseif (i >= 100) && (i <= 999)
                load([recsDir1(k,:),'/',cards(j,:),'/sample','00',num2str(i),'.mat'],'temp','rec');
            elseif (i >= 1000) && (i <= 9999)
                load([recsDir1(k,:),'/',cards(j,:),'/sample','0',num2str(i),'.mat'],'temp','rec');
            else
                load([recsDir1(k,:),'/',cards(j,:),'/sample',num2str(i),'.mat'],'temp','rec');
            end
    
            if (isempty(temp) == 0)
                temp_dat1(i) = temp(1);
            else
                temp_dat1(i) = 0.0;
            end
            
            noisy_mean1(i) = mean(rec(1:tSamp));
    
            noisy_std1(i) = std(rec(1:tSamp));
        
        end
    
        % ---- Getting Matched Filter Output ----
        load([testDir1(k,:),'/',cards(j,:),'_',cards(j,:),'.mat'],'m0');
        mfo_dat1 = m0.genOrig;

        % ---- Finding Available Records of Temperature Data ----
        temp_ind = find(temp_dat1~=0.0);
    
        temp_dat1 = temp_dat1(temp_ind);
        mfo_dat1 = mfo_dat1(temp_ind);
        noisy_mean1 = noisy_mean1(temp_ind);
        noisy_std1 = noisy_std1(temp_ind);
    
        % ---- Removing Outlier (MFO) ----
        [mfo_dat1,~] = outliers(mfo_dat1,5);
        nan_x = isnan(mfo_dat1);
    
        temp_dat1 = temp_dat1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);

        % ---- Removing Outlier (Temperature) ----
        [temp_dat1,~] = outliers(temp_dat1,5);
        nan_x = isnan(temp_dat1);
    
        temp_dat1 = temp_dat1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);
        
        % ---- Removing Outlier (Mean of Noisy Part) ----
        [noisy_mean1,~] = outliers(noisy_mean1,5);
        nan_x = isnan(noisy_mean1);
    
        temp_dat1 = temp_dat1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);
        
        % ---- Removing Outlier (STD of Noisy Part) ----
        [noisy_std1,~] = outliers(noisy_std1,5);
        nan_x = isnan(noisy_std1);
    
        temp_dat1 = temp_dat1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);
    
        % ---- Storing Data ----
        Temp_data = temp_dat1;
        MFO_data = mfo_dat1;
        Noisy_mean = noisy_mean1;
        Noisy_std = noisy_std1;
    
    save([ResDir,'/','type_1_s',num2str(k),'_',cards(j,:),'_data.mat'], ...
        'Temp_data','MFO_data','Noisy_mean','Noisy_std');

    end

end

%% ---- Getting Data of Dataset (Type 2) ----

disp('Getting Data for Dataset (Type 2):');

for k = 1:size(recsDir2,1)

    for j = 1:size(cards,1)

        % ---- Getting Temperature Data ----
        for i = 1:tRec2
    
            if mod(i,n2) == 0
                disp([cards(j,:) ' ' num2str(i/tRec2*100) '% complete (' num2str(i) '/' tRecStr2 ' recs analyzed)...']);
            end

            if (i >= 1) && (i <= 9)
                load([recsDir2(k,:),'/',cards(j,:),'/sample','0000',num2str(i),'.mat'],'temp','rec');
            elseif (i >= 10) && (i <= 99)
                load([recsDir2(k,:),'/',cards(j,:),'/sample','000',num2str(i),'.mat'],'temp','rec');
            elseif (i >= 100) && (i <= 999)
                load([recsDir2(k,:),'/',cards(j,:),'/sample','00',num2str(i),'.mat'],'temp','rec');
            elseif (i >= 1000) && (i <= 9999)
                load([recsDir2(k,:),'/',cards(j,:),'/sample','0',num2str(i),'.mat'],'temp','rec');
            else
                load([recsDir2(k,:),'/',cards(j,:),'/sample',num2str(i),'.mat'],'temp','rec');
            end
    
            if (isempty(temp) == 0)
                temp_dat2(i) = temp(1);
            else
                temp_dat2(i) = 0.0;
            end
            
            noisy_mean2(i) = mean(rec(1:tSamp));
    
            noisy_std2(i) = std(rec(1:tSamp));
        
        end
    
        % ---- Getting Matched Filter Output ----
        load([testDir2(k,:),'/',cards(j,:),'_',cards(j,:),'.mat'],'m0');
        mfo_dat2 = m0.genOrig;

        % ---- Finding Available Records of Temperature Data ----
        temp_ind = find(temp_dat2~=0.0);
    
        temp_dat2 = temp_dat2(temp_ind);
        mfo_dat2 = mfo_dat2(temp_ind);
        noisy_mean2 = noisy_mean2(temp_ind);
        noisy_std2 = noisy_std2(temp_ind);
    
        % ---- Removing Outlier (MFO) ----
        [mfo_dat2,~] = outliers(mfo_dat2,5);
        nan_x = isnan(mfo_dat2);
    
        temp_dat2 = temp_dat2(~nan_x);
        mfo_dat2 = mfo_dat2(~nan_x);
        noisy_mean2 = noisy_mean2(~nan_x);
        noisy_std2 = noisy_std2(~nan_x);

        % ---- Removing Outlier (Temperature) ----
        [temp_dat2,~] = outliers(temp_dat2,5);
        nan_x = isnan(temp_dat2);
    
        temp_dat2 = temp_dat2(~nan_x);
        mfo_dat2 = mfo_dat2(~nan_x);
        noisy_mean2 = noisy_mean2(~nan_x);
        noisy_std2 = noisy_std2(~nan_x);
        
        % ---- Removing Outlier (Mean of Noisy Part) ----
        [noisy_mean2,~] = outliers(noisy_mean2,5);
        nan_x = isnan(noisy_mean2);
    
        temp_dat2 = temp_dat2(~nan_x);
        mfo_dat2 = mfo_dat2(~nan_x);
        noisy_mean2 = noisy_mean2(~nan_x);
        noisy_std2 = noisy_std2(~nan_x);
        
        % ---- Removing Outlier (STD of Noisy Part) ----
        [noisy_std2,~] = outliers(noisy_std2,5);
        nan_x = isnan(noisy_std2);
    
        temp_dat2 = temp_dat2(~nan_x);
        mfo_dat2 = mfo_dat2(~nan_x);
        noisy_mean2 = noisy_mean2(~nan_x);
        noisy_std2 = noisy_std2(~nan_x);
    
        % ---- Storing Data ----
        Temp_data = temp_dat2;
        MFO_data = mfo_dat2;
        Noisy_mean = noisy_mean2;
        Noisy_std = noisy_std2;
    
    save([ResDir,'/','type_2_s',num2str(k),'_',cards(j,:),'_data.mat'], ...
        'Temp_data','MFO_data','Noisy_mean','Noisy_std');

    end

end

%% ---- Getting Data of Dataset (Type 3) ----

disp('Getting Data for Dataset (Type 3):');

for k = 1:size(recsDir3,1)

    for j = 1:size(cards,1)

        % ---- Getting Temperature and Voltage Data ----
        for i = 1:tRec1
    
            if mod(i,n1) == 0
                disp([cards(j,:) ' ' num2str(i/tRec1*100) '% complete (' num2str(i) '/' tRecStr1 ' recs analyzed)...']);
            end

            if (i >= 1) && (i <= 9)
                load([recsDir3(k,:),'/',cards(j,:),'/sample','0000',num2str(i),'.mat'],'volt_dat','temp','rec');
            elseif (i >= 10) && (i <= 99)
                load([recsDir3(k,:),'/',cards(j,:),'/sample','000',num2str(i),'.mat'],'volt_dat','temp','rec');
            elseif (i >= 100) && (i <= 999)
                load([recsDir3(k,:),'/',cards(j,:),'/sample','00',num2str(i),'.mat'],'volt_dat','temp','rec');
            elseif (i >= 1000) && (i <= 9999)
                load([recsDir3(k,:),'/',cards(j,:),'/sample','0',num2str(i),'.mat'],'volt_dat','temp','rec');
            else
                load([recsDir3(k,:),'/',cards(j,:),'/sample',num2str(i),'.mat'],'volt_dat','temp','rec');
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
            
            noisy_mean1(i) = mean(rec(1:tSamp));
    
            noisy_std1(i) = std(rec(1:tSamp));
        
        end
    
        % ---- Getting Matched Filter Output ----
        load([testDir3(k,:),'/',cards(j,:),'_',cards(j,:),'.mat'],'m0');
        mfo_dat1 = m0.genOrig;

        % ---- Finding Available Records of Temperature Data ----
        temp_ind = find(temp_dat1~=0.0);
    
        temp_dat1 = temp_dat1(temp_ind);
        sv_Tran1 = sv_Tran1(temp_ind);
        sv_IC1 = sv_IC1(temp_ind);
        mfo_dat1 = mfo_dat1(temp_ind);
        noisy_mean1 = noisy_mean1(temp_ind);
        noisy_std1 = noisy_std1(temp_ind);
    
        % ---- Removing Outlier (MFO) ----
        [mfo_dat1,~] = outliers(mfo_dat1,5);
        nan_x = isnan(mfo_dat1);
    
        temp_dat1 = temp_dat1(~nan_x);
        sv_Tran1 = sv_Tran1(~nan_x);
        sv_IC1 = sv_IC1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);

        % ---- Removing Outlier (Temperature) ----
        [temp_dat1,~] = outliers(temp_dat1,5);
        nan_x = isnan(temp_dat1);
    
        temp_dat1 = temp_dat1(~nan_x);
        sv_Tran1 = sv_Tran1(~nan_x);
        sv_IC1 = sv_IC1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);
    
        % ---- Removing Outlier (IC Supply Voltage) ----
        [sv_IC1,~] = outliers(sv_IC1,5);
        nan_x = isnan(sv_IC1);
    
        temp_dat1 = temp_dat1(~nan_x);
        sv_Tran1 = sv_Tran1(~nan_x);
        sv_IC1 = sv_IC1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);

        % ---- Removing Outlier (Transformer Supply Voltage) ----
        [sv_Tran1,~] = outliers(sv_Tran1,5);
        nan_x = isnan(sv_Tran1);
    
        temp_dat1 = temp_dat1(~nan_x);
        sv_Tran1 = sv_Tran1(~nan_x);
        sv_IC1 = sv_IC1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);
        
        % ---- Removing Outlier (Mean of Noisy Part) ----
        [noisy_mean1,~] = outliers(noisy_mean1,5);
        nan_x = isnan(noisy_mean1);
    
        temp_dat1 = temp_dat1(~nan_x);
        sv_Tran1 = sv_Tran1(~nan_x);
        sv_IC1 = sv_IC1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);
        
        % ---- Removing Outlier (STD of Noisy Part) ----
        [noisy_std1,~] = outliers(noisy_std1,5);
        nan_x = isnan(noisy_std1);
    
        temp_dat1 = temp_dat1(~nan_x);
        sv_Tran1 = sv_Tran1(~nan_x);
        sv_IC1 = sv_IC1(~nan_x);
        mfo_dat1 = mfo_dat1(~nan_x);
        noisy_mean1 = noisy_mean1(~nan_x);
        noisy_std1 = noisy_std1(~nan_x);
    
        % ---- Storing Data ----
        Temp_data = temp_dat1;
        Tran_volt = sv_Tran1;
        IC_volt = sv_IC1;
        MFO_data = mfo_dat1;
        Noisy_mean = noisy_mean1;
        Noisy_std = noisy_std1;
    
    save([ResDir,'/','type_3_s',num2str(k),'_',cards(j,:),'_data.mat'], ...
        'Temp_data','Tran_volt','IC_volt','MFO_data','Noisy_mean','Noisy_std');

    end

end
% -------------------------------