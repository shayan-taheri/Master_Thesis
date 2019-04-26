function [cm, aprs, cm_cnt] = build_fuse_cm_aprs()
% build a confusion matrix and aprs tables for the combined matched filter
% (ndss-journal)

% DATA INFORMATION
RecsLoc = '/Users/gerdes/dilon/';
RecsLoc_fn = [RecsLoc 'fn/'];
RecsLoc_fp = [RecsLoc 'fp/'];

% CARD INFORMATION
% cards to include in CM
cards = ['b4c1'; 'b4c2'; 'b4c3'];
% cards = [...
%     'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
%     'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; ...
%     ];
% cards = ['b6c1'; 'b6c2'; 'b6c3'];
% cards = ['b4c1 '; 'b4c2 '; 'b4c3 '; ...
%     'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
%     'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; ...
%     'b6c1 '; 'b6c2 '; 'b6c3 '];
% num of cards
n = length(cards(:,1));

% the confusion matrix
cm = zeros(n);
ns = zeros(n,1); %num of subject recs

% the aprs table (to build intra-model, use only same model in 'cards')
aprs = zeros(n,4);
cm_cnt = zeros(n); % cm of tn and fn counts

% build the cm...
for i = 1:n
    disp(['Filing in CM for ' deblank(cards(i,:)) '...']);
    
    % true negative...
    disp(['  Finding true negative rate...']);
    
    % load control
    load([RecsLoc_fp deblank(cards(i,:)) '_fp.mat'], 'fu_dfpi');
    
    % calc true negative rate using dfpi
    ns(i) = fu_dfpi.vars(2); %num of subject records
    cm_cnt(i,i) = sum(fu_dfpi.acpt) + fu_dfpi.vars(1);
    cm(i,i) = fu_dfpi.vars(3);
    
    % false negatives...
    for j = 1:n
        % no need for inter-card comparison (except for b6c2/3 against
        % b5...check by hand)
        if j ~= i && cards(i,2) == cards(j,2)
            disp(['  Finding false negative rate for ' deblank(cards(j,:)) '...']);
            % load subject records
            load([RecsLoc_fn deblank(cards(i,:)) '/' ...
                deblank(cards(j,:)) '_fn.mat'], 'fu_dfpi_fn');

            % calc false negative rate using thresholds from prediction 
            % intervals
            cm_cnt(i,j) = round(mean(fu_dfpi_fn.acpt_cnt));
            cm(i,j) = cm_cnt(i,j)/fu_dfpi_fn.vars(1);
        end
    end
end

% build aprs table...(we could redo this all as matrix operations)
ntrn = fu_dfpi.vars(1); %num of training recs (same for all cards)
for i = 1:n
    % find tn & fp...
    TN = cm_cnt(i,i);
    FP = ns(i) - ntrn - TN;
    
    % find fn & tp...
    FN = sum(cm_cnt(i,:)) - cm_cnt(i,i);
    TP = sum(ns) - ns(i) - FN;

    % accuracy
    aprs(i,1) = (TP+TN)/(TP+TN+FP+FN);
    % precision
    aprs(i,2) = TP/(TP+FP);
    % recall
    aprs(i,3) = TP/(TP+FN);
    % specifity
    aprs(i,4) = TN/(TN+FP);
end

% write cm to file as LaTeX table
cm2latex(cards, cm);
% write aprs table to file as LaTeX table
aprs2latex(cards, aprs);
end


function cm2latex(cards, cm)
disp('Writing CM to file with LaTeX notation (three decimal places)...');
% have to have leading zero
n = size(cm,1);

fid = fopen('cm_fuse.tab','w');
for i = 1:n
    card_c = deblank(cards(i,:));
    fprintf(fid, '%s & ', card_c);
    
    for j = 1:n
        if j == n
            cm_pv(fid,cm(i,j));
            fprintf(fid, '\\\\ \\hline\n');
        else
            cm_pv(fid,cm(i,j));
            fprintf(fid, '& ');
        end
    end 
end
fclose(fid);
end

% print val of cm
function cm_pv(fid, val)
val = round2(val,.001);

if val >= 0.001
    s = num2str(val,'%.3f');
    fprintf(fid, '%s ', s);
else
    fprintf(fid, '0 ');
end
end

function aprs2latex(cards, aprs)
disp('Writing aprs table to file with LaTeX notation (three decimal places)...');
% have to have leading zero
n = size(aprs,1);

fid = fopen('aprs_fuse.tab','w');
for i = 1:n
    card_s = deblank(cards(i,:));
    for j = 1:4
        if j == 1
            fprintf(fid, '%s & ', card_s);
            aprs_pv(fid,aprs(i,j));
            fprintf(fid, '& ');
        elseif j == 4
            aprs_pv(fid,aprs(i,j));
            fprintf(fid, '\\\\ \n');
        else
            aprs_pv(fid,aprs(i,j));
            fprintf(fid, '& ');
        end
    end
end
fclose(fid);
end

% print val of aprs
function aprs_pv(fid, val)
val = round2(val,.001);

s = num2str(val,'%.3f');
fprintf(fid, '%s ', s);
end