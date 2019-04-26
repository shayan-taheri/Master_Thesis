function g_fp(conCard,conDir, fpDir)
% function for the fp case: load data and evaluate it
% conCard: name of control card
% conDir: directory of control data
% fpDir: where to save fp results

% remember to execute matlabpool(n), where n is number of cores, at command
% window before running

% directory delimiter
DD = '/'; 

% load filter outputs
load([conDir DD conCard '_' conCard]);
m0con = m0;

% number of tests; each test was normed, hence the two factor
% bpf = 210; gen = 1; trimL = 12; trimU = 12; trimLRef = 12; trimURef = 12; fary = 10;
% nt = 2*(bpf + gen + trimL + trimU + trimLRef + trimURef + fary);
nt = 1;

% *** Important *** number of control records, nominal = num 0f records to load
nc = 1000; % nc = 10000;
% m0: all tests
m0 = zeros(nc,nt);
% num of outputs to use for traning in fp, ind tests
n = 25;
%num of outputs to calc thresholds for, ind tests
m = 20;

% put outputs in m0 matrix: all tests followed by their norm
m0(:,1:1) = m0con.gen';
% m0(:,2:13) = m0con.trimL';
% m0(:,14:25) = m0con.trimU';
% m0(:,26:37) = m0con.trimLRef';
% m0(:,38:49) = m0con.trimURef';
% m0(:,50:59) = m0con.fary';
% m0(:,60:269) = m0con.bpf(1:bpf,:)';
% m0(:,270:270) = m0con.genNorm';
% m0(:,271:282) = m0con.trimLNorm';
% m0(:,283:294) = m0con.trimUNorm';
% m0(:,295:306) = m0con.trimLRefNorm';
% m0(:,307:318) = m0con.trimURefNorm';
% m0(:,319:328) = m0con.faryNorm';
% m0(:,329:538) = m0con.bpfNorm(1:bpf,:)';

% use only 'effective' tests
% load(['eft/eft_' Recs '.mat']);
% m0 = m0(:,eft_c);
% m0 = m0(:,[1 328]); %b4
% m0 = m0(:,[1 162]); %b4
% m0 = m0(:,[1 168 477]); %b5
% m0 = m0(:,[1 50]); %b6
% nt = length(m0(1,:));

% delete old vars
clear m0con;

% remove bad recs
m0_orig = m0(:,1); %keep one orig output for debug
ind = badRec(m0(:,1));
nb = length(ind); %num of bad records
% remove bad recs from m0 completely (shrink m0)
m0(ind,:) = 0;
m0 = m0(m0(:,1)~=0,:);

% *** Important *** num of control records, actual
nc = 1000-nb; % nc = 10000-nb;
%num of m-period threshold sets, ind tests
th_p = ceil((nc-n)/m);
% thresholds: num th periods by upper/lower bounds by num tests
th = zeros(th_p,2,nt);
% indicates record status
%   0: rec failed test
%   1: rec passed test
acpts = zeros(nt,nc,'int8');

% run fp routine on individual tests: norm tests can have wider th than regular
% ones (see fp func)
tic;
disp('fp tests...');
parfor i = 1:nt
%    disp(['  test: ' num2str(i)]);
    fp_st = fp(m0(:,i),0);
    acpts(i,:) = fp_st.acpt';
    th(:,:,i) = fp_st.th;
end
toc;
% % norm of tests
% tic;
% disp('norm of tests...');
% parfor i = 212:421
% %     disp(['test: ' num2str(i)]);
%     fp_st = fp(m0(:,i),1);
%     acpts(i,:) = fp_st.acpt';
%     th(:,:,i) = fp_st.th;
% end
% toc;

% basic statistics on number of accepts, individual tests
a = sum(acpts(:,n+1:nc)); %acpts per rec
a_t = sum(acpts(:,n+1:nc),2); %acpts per test

% % fusion of data, prediction interval approach
% tic;
% disp('fusion...');
% fu_dfpi = fuse_fp_dfpi(a,acpts(1,n+1:nc));
% toc;
% % % fusion of data, clustering
% % fu_clu = fuse_fp_clu(acpts);

% % find which tests aren't effective: moved to g_eft3
% tind=1:nt;
% neft = tind(a_t/(nc-n) < .98); %which tests were ineffective; may have to increase this, to say .99 (better fusion?)
% neft_tnr = a_t(neft)/(nc-n); %how ineffective were 'ineffective' tests

% save almost everything
clear m0;
save([fpDir DD conCard '_' conCard]);
