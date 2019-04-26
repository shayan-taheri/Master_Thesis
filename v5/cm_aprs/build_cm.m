function cm = build_cm()
% build a confusion matrix for the generic matched filter (ndss-journal)

% DATA INFORMATION
RecsLoc = '/Users/rgerdes/Documents/dilon/research/data/results/';

% CARD INFORMATION
% cards to include in CM
cards = [ ...
    'b4c1 '; 'b4c2 '; 'b4c3 '; 'b4c4 '; 'b4c5 '; ...
    'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
    'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; ...
    'b6c1 '; 'b6c2 '; 'b6c3 '; 'b6c4 '];
% num of cards
n = length(cards(:,1));

% the confusion matrix
cm = zeros(n);

% build the cm...
for i = 1:n
    card_c = deblank(cards(i,:));
    disp(['Filing in CM for ' card_c '...']);
    
    % true negative...
    disp(['  Finding true negative rate...']);
    
    % load control records (generic filter output)
    load([RecsLoc 'm0' card_c '.mat'], 'm0ref');
    m0 = m0ref.gen;
    
    % scrub records
    ind = badRec(m0); %get indices of bad recs
    % remove bad recs from m0 completely (shrink m0)
    m0(ind) = 0;
    m0 = m0(m0~=0,:);

    % calc true negative rate using prediction intervals
    fp_st = fp(m0,1);
    cm(i,i) = sum(fp_st.acpt)/(fp_st.vars(4) - fp_st.vars(1));

    % false negatives...
    for j = 1:n
        if j ~= i
            card_t = deblank(cards(j,:));
            disp(['  Finding false negative rate for ' card_t '...']);
            % load subject records (generic filter output)
            Recs_int = ['m0int_' card_t];
            load([RecsLoc 'm0' card_c '.mat'], Recs_int);
            eval(['m0int = ' Recs_int '.gen;']);

            % scrub records
            ind = badRec(m0int); %get indices of bad recs
            % remove bad recs from m0int completely (shrink m0int)
            m0int(ind) = 0;
            m0int = m0int(m0int ~= 0,:);

            % calc false negative rate using thresholds from prediction 
            % intervals
            fn_st = fn(m0int,fp_st.th);
            cm(i,j) = mean(fn_st.acpt_cnt)/fn_st.vars(1);
%             % debug
%             figure;plot(m0);hold on;plot(m0int,'r');
%             title([card_c ' vs. ' card_t]);
        end
    end
end

% write cm to file as LaTeX table
mat2latex(cards,cm);

function mat2latex(cards,cm)
disp('Writing CM to file using LaTeX notation (three decimal places)...');
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

% print val of cm
function cm_pv(fid, val)
val = round2(val,.001);

if val >= 0.001
    s = num2str(val,'%.3f');
    fprintf(fid, '%s ', s);
else
    fprintf(fid, '0 ');
end