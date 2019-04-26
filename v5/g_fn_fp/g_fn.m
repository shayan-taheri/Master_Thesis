function g_fn(conCard, subCard, subDir, fnDir, th)
% wrapper for the fn case: load data and evaluate it
% conCard: name of control card
% subCard: name of subject card
% subDir: directory of subject data
% fnDir: where to save fn results
% th: thresholds of individual tests (num th periods by upper/lower bounds
%     by num tests)

% remember to execute matlabpool(n), where n is number of cores, at command
% window before running

% directory delimiter
DD = '/'; 

% load filter outputs
load([subDir DD conCard '_' subCard]);
m0sub = m0;

% number of tests; each test was normed, hence the two factor
% bpf = 210; gen = 1; trimL = 12; trimU = 12; trimLRef = 12; trimURef = 12; fary = 10;
% nt = 2*(bpf + gen + trimL + trimU + trimLRef + trimURef + fary);
nt = 1;

% *** Important *** number of subject outputs, nominal = num 0f records to load
ns = 1000; % ns = 10000;
% m0: all tests
m0 = zeros(ns,nt);

% put outputs in m0 matrix: all tests followed by their norm
m0(:,1:1) = m0sub.gen';
% m0(:,2:13) = m0sub.trimL';
% m0(:,14:25) = m0sub.trimU';
% m0(:,26:37) = m0sub.trimLRef';
% m0(:,38:49) = m0sub.trimURef';
% m0(:,50:59) = m0sub.fary';
% m0(:,60:269) = m0sub.bpf(1:bpf,:)';
% m0(:,270:270) = m0sub.genNorm';
% m0(:,271:282) = m0sub.trimLNorm';
% m0(:,283:294) = m0sub.trimUNorm';
% m0(:,295:306) = m0sub.trimLRefNorm';
% m0(:,307:318) = m0sub.trimURefNorm';
% m0(:,319:328) = m0sub.faryNorm';
% m0(:,329:538) = m0sub.bpfNorm(1:bpf,:)';


% % use only 'effective' tests
% % load(['eft/eft_' Recs '.mat']);
% % m0 = m0(:,eft_c);
% % nt = length(m0(1,:));
% % m0 = m0(:,[1 328]); %b4
% % m0 = m0(:,[1 168 477]); %b5
% m0 = m0(:,[1 50]); %b6
% nt = length(m0(1,:));

% delete old vars
clear m0sub;
% remove bad recs
m0_orig = m0(:,1); %keep one orig output for debug
ind = badRec(m0(:,1));
nb = length(ind); %num of bad records
% remove bad recs from m0 completely (shrink m0)
m0(ind,:) = 0;
m0 = m0(m0(:,1)~=0,:);

% *** Important *** num of subject records, actual
ns = 1000-nb; % ns = 10000-nb;
%num of m-period threshold sets, ind tests, for control
th_p = length(th(:,1,1));
% acpt_cnt: indicates num subject rec within thresholds calc for control
%   0: no records within th
%   n: n-rec(s) within th; the average of this is our FN count, for ind
%      tests
acpt_cnt = zeros(th_p,nt,'uint16');

% % debug of fn_s()
% fn_s_cnt = zeros(th_p,nt,'uint16');
% % fn_st.acpt_cnt(i) = sum(fn_st.acpts(:,i));
% tic;
% disp('fn_s: all tests...');
% parfor i = 1:th_p
%     fn_s_st = fn_s(m0,th(i,:,:));
%     fn_s_cnt(i,:) = sum(fn_s_st.acpts);
% end
% toc;    

% run fn routine on individual tests
tic;
disp('fn: individual tests...');
parfor i = 1:nt
%    disp(['  test: ' num2str(i)]);
    fn_st = fn(m0(:,i), th(:,:,i));
    acpt_cnt(:,i) = fn_st.acpt_cnt;
end
toc;

% basic statistics on number of accepts, individual tests
% fnc_s = mean(fn_s_cnt);
fnc = mean(acpt_cnt); %false-negative count, per test
% % debug
% figure;plot(fnc)
% hold on;
% plot(fnc_s,'rx');
% fnc(1:10)
% fnc_s(1:10)
% sum(fnc_s == fnc)

% % fusion of data, prediction interval approach
% tic;
% disp('fn: dfpi fusion');
% fu_dfpi_fn = fuse_fn_dfpi(m0,th,fth);
% toc;
% % % fusion of data, clustering
% % tic;
% % disp('fn: clu fusion');
% % fu_clu_fn = fuse_fn_clu(m0,th,fu_clu.th,fu_clu.idx,fu_clu.p,fu_clu.P);
% % toc;

% disp(['DFPI FNC: ' num2str(mean(fu_dfpi_fn.acpt_cnt))]);
% % disp(['CLU FNC: ' num2str(mean(fu_clu_fn.acpt_cnt))]);
% % % debug
% % figure;plot(fu_dfpi_fn.acpt_cnt);
% % figure;plot(fu_clu_fn.acpt_cnt);

% % find effective tests: moved to g_eft3
% tind=1:nt;
% eft = tind(mean(acpt_cnt)/ns < .50); %which tests were effective; may have to increase this
% eft_tpr = mean(acpt_cnt(:,eft))/ns; %how effective were 'effective' tests

% save nearly everything
clear m0;
save([fnDir DD conCard '_' subCard]);
