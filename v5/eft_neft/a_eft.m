% Analyse effective tests for the fn, all cards.  Effective tests 
% cutoff is specified herein (fn).  Requires access to fn 
% calculations for ALL tests.  Output is histogram of effective tests 

% DATA INFORMATION
% RecsLoc = '/Users/gerdes/dilon/';
RecsLoc = '/Volumes/mms0/dilon/';
RecsLoc_fn = [RecsLoc 'fn/'];

% CARD INFORMATION
% cards to include in CM
% cards = ['b4c1'; 'b4c2'; 'b4c3'];
cards = ['b4c1 '; 'b4c2 '; 'b4c3 '; ...
    'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
    'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; ...
    'b6c1 '; 'b6c2 '; 'b6c3 '];
% num of cards
n = length(cards(:,1));

% test index
tind=1:526;
% num rec used for training (fp)
ntrn = 25;
% occurence of eft tests
eft_cnt = zeros(size(tind,2),1); 

% the acceptable false-negative rate (lower means more strict)
fnr = .1;
% effective tests (against all subjects)
eft_s = [0 1]; %we include the first test because it is excluded in fn case (see below)


% find effective tests for fn case...
for i = 1:n
    card_c = deblank(cards(i,:));
    for j = 1:n
        % load effective tests for subjects
        card_s = deblank(cards(j,:));
        % no need for inter-card comparison (except for b6c2/3 against
        % b5...check by hand)
        if i ~= j && card_c(2) == card_s(2)
            load([RecsLoc_fn card_c '/' card_s '_fn.mat'],'acpt_cnt','ns');
            eft = tind(mean(acpt_cnt)/ns < fnr); %which tests were effective; may have to decrease this
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

save('eft_fn.mat','eft_s','eft_cnt');

% % find how often the effective tests are effective (count the num of cards
% % a particular test is effective for)
% 
% % num of effective tests
% eft_n = size(eft_s,2);
% % occurence of eft tests
% eft_cnt = zeros(eft_n,1); 
% 
% % count eft occurences...
% for i = 1:n
%     card_c = deblank(cards(i,:));
%     for j = 1:n
%         % load effective tests for subjects
%         card_s = deblank(cards(j,:));
%         if i ~= j && card_c(2) == card_s(2)
%             load([RecsLoc_fn card_c '/' card_s '_fn.mat'],'acpt_cnt','ns');
%             eft = (mean(acpt_cnt)/ns < fnr); %which tests were effective
%             
%             if size(eft,2) > 0 && eft(1) ~= 1 %if gen mf effective, ignore all others
%                 eft_cnt(eft) = eft_cnt(eft) + 1;
%             end
%         end
%     end
% end