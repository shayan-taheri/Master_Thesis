% *** Storing MFO and Noisy Part using Extracted Data ***
%% ---- Preparing Environment ----
clear all;
close all;

cards = char('b5c1','b5c2','b5c3','b5c5','b5c7');

loadDir = '/media/SHAYAN_HDD/Results/Analysis_9';

resDir = '/media/SHAYAN_HDD/Results/Analysis_10';

% ---- Loading Extracted Data ----

for i = 1:3 % Iteration in Types
    
    for j = 1:2 % Iteration in Datasets if each Type
        
        set_index = 2*(i-1)+j;
        
        for k = 1:size(cards,1) % Iteration in Cards of each Dataset
            
            load([loadDir,'/type_',num2str(i),'_s',num2str(j),'_', ...
                cards(k,:),'_data.mat']);
            
            eval(['dataset',num2str(set_index),'.',cards(k,:),'.','MFO_data = MFO_data']);
            eval(['dataset',num2str(set_index),'.',cards(k,:),'.','Noisy_mean = Noisy_mean']);
            eval(['dataset',num2str(set_index),'.',cards(k,:),'.','Noisy_std = Noisy_std']);
            
        end
        
    end
    
end

% ---- Storing Desired Data ----
save([resDir,'/','mfo_noisy_data.mat'], ...
        'dataset1','dataset2','dataset3','dataset4','dataset5','dataset6');
% --------------------------------------------------------------