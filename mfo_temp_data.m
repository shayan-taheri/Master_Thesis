% *** Storing MFO and Temperature Data in ".mat" File ***
% ---- Preparing Environment ----
clear all;
close all;
card = 'b5c2';
testDir = '/media/SHAYAN_HDD/Results/Collection_4/test_res_reg/1';
recsDir = ['/media/SHAYAN_HDD/Data/collection_4/',card];
ResDir = '/media/SHAYAN_HDD/Results/Analysis_7/mfo_temp_data';
tRec = 5308; % Number of Records
% ---- Fetching Temperature Data ----
temp_orig = zeros(1,tRec); % Vector for Temperature Data
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
        temp_orig(i) = temp(1);
    else
        temp_orig(i) = NaN;
    end
    
end

% Finding Good Records: The ones that have temperature data.
nan_ind = isnan(temp_orig);
temp_trim = temp_orig(~nan_ind);

% ---- Matched Filter Output ----
mfo_orig = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
mfo_orig = m0.genOrig;
mfo_trim = mfo_orig(~nan_ind); % Correspondence of Temperature Data

% ---- Storing Data ----
save([ResDir,'/',card,'_data.mat'],'mfo_orig','temp_orig','mfo_trim','temp_trim');
% -------------------------------