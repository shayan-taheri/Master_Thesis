function [cm, aprs, cm_cnt] = build_cm_aprs()
% build a confusion matrix and aprs tables for the generic matched filter
% (ndss-journal)

% DATA INFORMATION
RecsLoc = '/Users/rgerdes/Documents/dilon/research/data/results/';

% CARD INFORMATION
% cards to include in CM
cards = ['b4c1'; 'b4c2'; 'b4c3'; 'b4c4'; 'b4c5'];
% cards = [...
%     'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
%     'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; ...
%     ];
% cards = ['b6c1'; 'b6c2'; 'b6c3'; 'b6c4'];
cards = [...
    'b4c1 '; 'b4c2 '; 'b4c3 '; 'b4c4 '; 'b4c5 '; 'b4c6 ';...
    'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
    'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; 'b5c11';...
    'b6c1 '; 'b6c2 '; 'b6c3 '; 'b6c4 '; 'b6c5 '];
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
    card_c = deblank(cards(i,:));
    disp(['Filing in CM for ' card_c '...']);
    
    % true negative...
    disp(['  Finding true negative rate...']);
    
    % load control records (generic filter output)
    load([RecsLoc '/m0' card_c '.mat'], 'm0ref');
    m0 = m0ref.gen;
    
    % scrub records
    ind = badRec(m0); %get indices of bad recs
    % remove bad recs from m0 completely (shrink m0)
    m0(ind) = 0;
    m0 = m0(m0~=0,:);

    % calc true negative rate using prediction intervals
    if card_c(2) ~= '6'
        fp_st = fp(m0,0); %b4/5
    else
        fp_st = fp(m0,1); %b6
    end
    ns(i) = fp_st.vars(4);
    cm_cnt(i,i) = sum(fp_st.acpt(fp_st.vars(1)+1:end));
    cm(i,i) = cm_cnt(i,i)/(ns(i) - fp_st.vars(1));
    
    % false negatives...
    for j = 1:n
        if j ~= i
            disp(['  Finding false negative rate for ' deblank(cards(j,:)) '...']);
            % load subject records (generic filter output)
            Recs_int = ['m0int_' deblank(cards(j,:))];
            load([RecsLoc '/m0' deblank(cards(i,:)) '.mat'], Recs_int);
            eval(['m0int = ' Recs_int '.gen;']);

            % scrub records (really, it would be more efficient to do this
            % just once for all records; i.e., load them and scrub them 
            % just once)
            ind = badRec(m0int); %get indices of bad recs
            % remove bad recs from m0int completely (shrink m0int)
            m0int(ind) = 0;
            m0int = m0int(m0int ~= 0,:);

            % calc false negative rate using thresholds from prediction 
            % intervals
            fn_st = fn(m0int,fp_st.th);
            cm_cnt(i,j) = round(mean(fn_st.acpt_cnt));
            cm(i,j) = cm_cnt(i,j)/fn_st.vars(1);
%             % debug
%             figure;plot(m0);hold on;plot(m0int,'r');
%             title([deblank(cards(i,:)) ' vs. ' deblank(cards(j,:))]);
        end
    end
end

% build aprs table...(we could redo this all as matrix operations)
ntrn = fp_st.vars(1); %num of training recs (same for all cards)
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

fid = fopen('cm.tab','w');
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

fid = fopen('aprs.tab','w');
for i = 1:n
    card_c = deblank(cards(i,:));
    fprintf(fid, '%s & ', card_c);
    
    for j = 1:4
        if j == 4
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