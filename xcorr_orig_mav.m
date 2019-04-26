% *** Performing Cross Correlation on Original Data ***
%% **** Preparing Environment and Getting Temperature Data ****
clear all;
close all;
card = 'b5c7';
recsDir = '/media/SHAYAN_HDD/Data/collection_2/b5c7';
testDir = '/media/SHAYAN_HDD/Results/Collection_2/test_res_mean/name';
tRec = 10000; % Number of Records
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

%% **** Matched Filter Output ****
y = zeros(1,tRec); % Vector for Matched Filter Output
load([testDir,'/',card,'_',card,'.mat'],'m0');
y = m0.genOrig;
y = y(g_ix); % Using matched filter outputs of good records.

% ---- Performing Cross Correlation ----
[Outputs,Lags] = xcorr(x,y);
[~,Index] = max(abs(Outputs));

fig_id = figure();
plot(Lags,Outputs);
xlabel('Lag', 'FontSize', 14);
ylabel('Cross Correlation Output', 'FontSize', 14);
title('Cross Correlation Analysis for MFO', 'FontSize', 16);
set(gca, 'fontsize', 12);
fig_leg = legend (['Maximum Output at Lag: ',num2str(Lags(Index))]);
set(fig_leg,'FontSize',8);
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_2/cross_res/',card,'_cross_orig','.pdf']);

% ---- Performing Sample Cross Correlation Function ----
fig_id = figure();
crosscorr(x,y);
xlabel('Lag', 'FontSize', 14);
ylabel('Sample Cross Correlation', 'FontSize', 14);
title('Sample Cross Correlation Function for MFO', 'FontSize', 16);
set(gca, 'fontsize', 12);
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_2/cross_res/',card,'_sample_orig','.pdf']);

%% **** Mean of Absolute Value (MAV) ****
z = zeros(1,tRec); % Vector for Mean of Absolute Value
z = m0.recs_AbsMean;
z = z(g_ix); % Using MAV of good records.

% ---- Performing Cross Correlation ----
[Outputs,Lags] = xcorr(x,z);
[~,Index] = max(abs(Outputs));

fig_id = figure();
plot(Lags,Outputs);
xlabel('Lag', 'FontSize', 14);
ylabel('Cross Correlation Output', 'FontSize', 14);
title('Cross Correlation Analysis for MAV', 'FontSize', 16);
set(gca, 'fontsize', 12);
fig_leg = legend (['Maximum Output at Lag: ',num2str(Lags(Index))]);
set(fig_leg,'FontSize',8);
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_2/cross_res/',card,'_cross_mav','.pdf']);

% ---- Performing Sample Cross Correlation Function ----
fig_id = figure();
crosscorr(x,z);
xlabel('Lag', 'FontSize', 14);
ylabel('Sample Cross Correlation', 'FontSize', 14);
title('Sample Cross Correlation Function for MAV', 'FontSize', 16);
set(gca, 'fontsize', 12);
saveas(fig_id,['/media/SHAYAN_HDD/Results/Collection_2/cross_res/',card,'_sample_mav','.pdf']);

% -------------------------------------------------------