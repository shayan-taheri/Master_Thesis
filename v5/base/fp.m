function fp_st = fp(m0,norm)
% fp(): determine status and threshold of records
% m0: filter output, for a single test
% norm: 1: recs are norm'd, 0: recs aren't norm'd
% fp_st: the false positive struct, see below

% variables
n = 25; %num of outputs to use for traning
m = 20; %num of outputs to calc thresholds for
if norm == 1 %this option is no longer about norm, but used for b6 exclusively
    r = 3.084;%range parameter for n and m, with alpha = .10; b6
else
%     r = 3.084;%range parameter for n and m, with alpha = .10; b6
    r = 3.397;%3.084;%range parameter for n and m, with alpha = .05; b4/5
end
nc = length(m0); %num of control records
th_p = ceil((nc-n)/m); %num of m-period threshold sets

% define false positives struct
% th: thresholds for each m records
% acpt: indicates record status
%   0: outside th
%   1: within th
%   -1: training data (don't include in fp calculations)
fp_st = struct(...
    'vars', [n,m,r,nc,th_p,norm], ...
    'th', zeros(th_p,2), ...
    'acpt', zeros(nc,1,'int8'));

% assign first n records to training vector
trn = m0(1:n);
% get first threshold set
fp_st.th(1,:) = g_th(trn,r);
% first n records have unique status
fp_st.acpt(1:n) = -1;
% determine thresholds and status for remaining records
for i = 1:th_p-1
    % rec indices
    ind_b = n+(i-1)*m+1; ind_e = n+i*m;
    % rec to test
    tst = m0(ind_b:ind_e);
    % rec within th
    acpt = (fp_st.th(i,1) <= tst) & (tst <= fp_st.th(i,2));
    % store rec status
    fp_st.acpt(ind_b:ind_e) = acpt;
    % get new training records and thresholds
    trn = g_trn(trn,tst,acpt,n,m);
    fp_st.th(i+1,:) = g_th(trn,r);
end
% will be nc-n+(th_p-1)*m rec left, det status;
if ind_e < nc
    tst = m0(ind_e+1:nc);
    fp_st.acpt(ind_e+1:nc) = ...
        (fp_st.th(i+1,1) <= tst) & (tst <= fp_st.th(i+1,2));
end

function th = g_th(trn,r)
th = zeros(1,2);
% if we had bad outputs, we should use nanmean/nanstd here...we should have
% removed those from the input though; this way we get an error
u_trn = mean(trn);
s_trn = std(trn);

% % evfit...
% p = evfit(trn);
% [u_trn,v_trn] = evstat(p(1),p(2));
% s_trn = sqrt(v_trn);

th(1) = u_trn - r*s_trn;
th(2) = u_trn + r*s_trn;

function trn = g_trn(trn,tst,acpt,n,m)
% old trn rec to retain
rtn = n - m + sum(acpt==0);
% use most recent outputs of old trn data for new trn data
trn(1:rtn) = trn(n-rtn+1:n);
if rtn < n %ensure that new recs were accepted...    
    % indices of acpt rec
    ind = find(acpt~=0);
    trn(rtn+1:n) = tst(ind);
end
