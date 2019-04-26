function [cm, aprs, cm_cnt] = g_cm_aprs(tst)
% build a confusion matrix and aprs tables for arbitary test; FP and FN
% values for test must alread by calculated; based on build_cm_aprs
% (which calcs fp & fn on the fly)
% tst: test to calc cm & aprs values for

% LOCATION of input/output data
% dir for results of fp tests
fpDir = '/home/shayan/pli/data/SP2/fn_fp/fp_test';
% dir for results of fn tests
fnDir = '/home/shayan/pli/data/SP2/fn_fp/fn_test';

% directory delimiter
DD = '/'; 

% CARD INFORMATION
% cards to include in CM
% cards = char('b4c1', 'b4c2', 'b4c3', 'b4c4', 'b4c5', 'b4c6');
% cards = [...
%     'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
%     'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; ...
%     ];
% cards = ['b6c1'; 'b6c2'; 'b6c3'; 'b6c4'];
cards = char(...
    'b4c1','b4c2','b4c3','b4c4','b4c5','b4c6',...
    'b5c1','b5c2','b5c3','b5c4','b5c5', ...
    'b5c6','b5c7','b5c8','b5c9','b5c10', ...
    'b6c1','b6c2','b6c3','b6c4','b6c5', ...
    'b6c6','b6c7','b6c8','b6c9','b6c10');
% num of cards
n = size(cards,1);

% the confusion matrix
cm = zeros(n); %rates, not counts
ns = zeros(n,1); %num of subject recs

% the aprs table (to build intra-model, use only same model in 'cards')
aprs = zeros(n,4);
cm_cnt = zeros(n); % cm of tn and fn counts

% build the cm(s)...
for i = 1:n
    conCard = deblank(cards(i,:));
    disp(['filing in CM for ' conCard '...']);
    % load fp data
    fp = load([fpDir DD conCard '_' conCard],'a_t','nc','n','nb');
    
    % true negative...
    disp(['  finding true negative rate...']);
    ns(i) = fp.nc;
    cm_cnt(i,i) = fp.a_t(tst);
    cm(i,i) = cm_cnt(i,i)/(ns(i) - fp.n);
    
    % false negatives...
    disp(['  finding false negative rates...']);
    for j = 1:n
        if j ~= i
            subCard = deblank(cards(j,:));
%             disp(['  finding false negative rate for ' conCard ' vs. ' subCard '...']);
            % load fn data
            fn = load([fnDir DD conCard '_' subCard],'acpt_cnt','ns');
            
            % false negative...
            % some cards may have ns = 0, which means that all recs were bad;
            % really, this means no overlap; keep cm at zero for them (no
            % division by zero)
            if fn.ns ~= 0
                cm_cnt(i,j) = round(mean(fn.acpt_cnt(:,tst)));
                cm(i,j) = cm_cnt(i,j)/fn.ns;
            end
            
        end
    end
end

% build aprs table...(we could redo this all as matrix operations)
ntrn = fp.n; %num of training recs (same for all cards)
% determine indices of cm that correspond to each model (mf and ml)
[b,mf,x] = unique(cards(:,1:2),'rows','first');
[b,ml,x] = unique(cards(:,1:2),'rows','last');
for i = 1:n
    % want intra-model calcs only; set limits on which parts of cm are used
    % depending upon model...the following does this by determining the
    % model number that corresponds to device i and then finding where that
    % model begins and ends in the cm
    m = cards(i,1:2)

    for j = 1:size(b,1) %b is the unique list of models
        if m == b(j,1:2) %if model(i) equals list of unique models, use indices for cm computed above
            cm_f = mf(j)
            cm_l = ml(j)
        end
    end
    
    % find tn & fp...
    TN = cm_cnt(i,i);
    FP = ns(i) - ntrn - TN;
    
    % find fn & tp...
    FN = sum(cm_cnt(i,cm_f:cm_l)) - cm_cnt(i,i);
    TP = sum(ns(cm_f:cm_l)) - ns(i) - FN;

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

fid = fopen('/home/shayan/pli/data/SP2/fn_fp/cm.tab','w');
% write table header
cm_header(fid,cards);
% print subj values
for i = 1:n
    conCard = deblank(cards(i,:)); 
    conCard(1) = 'm'; %move away from bXcY to mXcY
    fprintf(fid, '%s & ', conCard);
    
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

% write cm table header for latex
function cm_header(fid,cards)
cards(:,1) = 'm';
n = size(cards,1); %num of cards
n_str = num2str(n);
m = unique(cards(:,1:2),'rows'); %card models

fprintf(fid,['\\begin{tabular}{|l|' repmat('c|',1,n) '}\n']);
fprintf(fid,'\\hline\n');
fprintf(fid,['& \\multicolumn{' n_str '}{c}{Subject} \\vline \\\\\n']);
for i = 1:size(m,1)
    nm = length(strmatch(m(i,:),cards)); %num of cards in model class
    fprintf(fid,['& \\multicolumn{' num2str(nm) '}{c}{' m(i,:) '} \\vline ']);
end
fprintf(fid,'\\\\\n');
fprintf(fid,'Control');
for i = 1:n
    fprintf(fid,[' & ' deblank(cards(i,3:end))]);
end
fprintf(fid,'\\\\ \\hline\n');
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

fid = fopen('/home/shayan/pli/data/SP2/fn_fp/aprs.tab','w');
for i = 1:n
    conCard = deblank(cards(i,:));
    conCard(1) = 'm';
    fprintf(fid, '%s & ', conCard);
    
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
% print means
fprintf(fid,'\\hline\n');
fprintf(fid,'\\textbf{mean} ');
for i = 1:4
    fprintf(fid,'& \\textbf{');
    aprs_pv(fid,mean(aprs(:,i)));
    fprintf(fid,'} ');
end
fprintf(fid,'\\\\\n');

fclose(fid);
end

% print val of aprs
function aprs_pv(fid, val)
val = round2(val,.001);

s = num2str(val,'%.3f');
fprintf(fid, '%s ', s);
end