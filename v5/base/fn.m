function fn_st = fn(m0,th)
% fn(): determine status and threshold of test records, for subject
% m0: filter output, for a single test
% fp_vars: variables used in fp test
% th: thresholds calc from fp test
% fn_st: the false negative struct, see below

% variables
ns = length(m0); %num of subject records
th_p = length(th(:,1,1)); %num of m-period threshold sets, for control

% define false negative struct
% acpts: indicates record status
%   0: outside th
%   1: inside th
% acpt_cnt: indicates num subject rec within thresholds calc for control
%   0: no records within th
%   n: n-rec(s) within th; the average of this is our FN count
fn_st = struct(...
    'vars', [ns, th_p], ...
    'acpts', zeros(ns,th_p,'uint8'), ...
    'acpt_cnt', zeros(th_p,1,'uint16'));

% determine thresholds and status for subject records (should use matrix
% operations for this...)
for i = 1:th_p
    % rec within th
    fn_st.acpts(:,i) = ((th(i,1) <= m0) & (m0 <= th(i,2)));
    fn_st.acpt_cnt(i) = sum(fn_st.acpts(:,i));
end