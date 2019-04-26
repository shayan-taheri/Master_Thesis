% Get effective tests for the fn and fp case, all cards.  Effective tests 
% cutoff is specified herein (both fp and fn), with the added provision 
% that we don't include tests that are only effective for two or fewer
% cards.  Requires access to fn/fp calculations for ALL tests.  Output is 
% effective tests for each control card

% DATA INFORMATION
% dir for results of fn tests
fnDir = '/Users/rgerdes/Documents/dilon/research/eval_2_2/v5/results/2010-SPRING-01/1';
% dir for results of fp tests
fpDir = '/Users/rgerdes/Documents/dilon/research/eval_2_2/v5/results/2010-SPRING-01/1';
% directory delimiter
DD = '/';

% CARD INFORMATION
% cards = char('b4c1', 'b4c2', 'b4c3', 'b4c4', 'b4c5', 'b4c6');
% cards = char(...
%     'b5c1','b5c2','b5c3','b5c4','b5c5', ...
%     'b5c6','b5c7','b5c8','b5c9','b5c10','b5c11');
cards = char(...
    'b6c1','b6c2','b6c3','b6c4','b6c5', ...
    'b6c6','b6c7','b6c8','b6c9','b6c10');
% cards = char(...
%     'b4c1','b4c2','b4c3','b4c4','b4c5','b4c6',...
%     'b5c1','b5c2','b5c3','b5c4','b5c5', ...
%     'b5c6','b5c7','b5c8','b5c9','b5c10','b5c11',...
%     'b6c1','b6c2','b6c3','b6c4','b6c5', ...
%     'b6c6','b6c7','b6c8','b6c9','b6c10');
n = length(cards(:,1));

% test index
tind=1:538;
% num rec used for training (fp)
ntrn = 25;

% the acceptable false-negative rate (lower means more strict)
fnr = 0.5;
% the acceptable true-negative rate (higher is stricter)
tnr = 0.90;
% effective tests (against all subjects)
eft_s = [0 1]; %we include the first test because it is excluded in fn case (see below)
% occurence of eft_s tests
eft_cnt = zeros(size(tind,2),1);
eft_cnt(1) = 10; %gen mf is effective by definition
% min num of cards test must be effective for us to use it
% eft_min = 4; %b4
% eft_min = 6; %b5
eft_min = 6; %b6
% ineffective tests (for the control)
neft_c = 0;


% fn case...
for i = 1:n
    conCard = deblank(cards(i,:));
    for j = 1:n
        % load result of tests for subjects
        subCard = deblank(cards(j,:));
        % no need for inter-card comparison
        if i ~= j && conCard(2) == subCard(2)
            % load fn data as struct (don't want to overwrite existing
            % vars)
            fn = load([fnDir DD conCard '_' subCard],'acpt_cnt','ns');
            %which tests were effective; may have to decrease this
            eft = tind(mean(fn.acpt_cnt)/fn.ns < fnr);
            % if gen mf is effective against control, we don't include the
            % other tests
            if size(eft,2) > 0 && eft(1) ~= 1
                eft_s = unique([eft_s eft]);
                eft_cnt(eft) = eft_cnt(eft) + 1;
            end
        end
    end
end

disp('effective tests (subjects)...');
eft_s = setdiff(eft_s,0);
disp(eft_s);
disp(['  total: ' num2str(length(eft_s))]);
disp(' ');

% get rid of tests that are only effective for less than 'eft_min' cards
disp('effective tests (subjects) pruned...');
eft_s_p = tind(eft_cnt < eft_min);
eft_s = setdiff(eft_s,eft_s_p);
disp(eft_s);
disp(['  total: ' num2str(length(eft_s))]);
disp(' ');

% fp case...
for i = 1:n
    % load ineffective tests for control
    conCard = deblank(cards(i,:));
    % should load acpts as for fn case here
    fp = load([fpDir DD conCard '_' conCard],'a_t','nc','n','nb');
    % determine which tests were not effective for control    
    neft = tind(fp.a_t/(fp.nc - fp.n) < tnr);
    
    % effective tests = effective tests (fn) - ineffective tests (fp)
    eft_c = setdiff(eft_s, neft);
    
    disp(['total effective tests for ' conCard ': ' num2str(length(eft_c))]);
%     save(['eft/eft_' conCard '.mat'],'eft_c');
end