function fn_st = fn_s(m0,th)
% fn-s(): determine status of subject recs for a single th_p
% m0: filter output, all tests (ns-by-nt)
% th: thresholds calc from fp for all tests (1-by-2-by-nt)
% fn_st: the false negative struct, see below

% variables
ns = length(m0(:,1)); %num of subject records
nt = length(m0(1,:)); %num of tests

% define false negative struct
% acpts: indicates record status
%   0: outside th
%   1: inside th
fn_st = struct(...
    'vars', [ns,nt], ...
    'acpts', zeros(ns,nt,'uint8'));

% lower and upper th for period, in matrix form so we don't have to use a
% for loop when determining status for subject records
th = reshape(th,2,nt);
% repmat is slow, so we use bsxfun to accomplish this:
%   th_l = repmat(th(1,:),ns,1);
%   th_u = repmat(th(2,:),ns,1);
%   fn_st.acpts = ((th_l <= m0) & (m0 <= th_u));
th_l = bsxfun(@ge,m0,th(1,:));
th_u = bsxfun(@le,m0,th(2,:));

% rec within th
fn_st.acpts = (th_l & th_u);