% display overlap between effecitve and ineffective tests, all cards

% DATA INFORMATION
RecsLoc = '/Users/gerdes/dilon/';
RecsLoc_fn = [RecsLoc 'fn/'];
RecsLoc_fp = [RecsLoc 'fp/'];

% CARD INFORMATION
% cards to include in CM
% cards = ['b4c1'; 'b4c2'; 'b4c3'];
cards = ['b4c1 '; 'b4c2 '; 'b4c3 '; ...
    'b5c1 '; 'b5c2 '; 'b5c3 '; 'b5c4 '; 'b5c5 '; ...
    'b5c6 '; 'b5c7 '; 'b5c8 '; 'b5c9 '; 'b5c10'; ...
    'b6c1 '; 'b6c2 '; 'b6c3 '];
% num of cards
n = length(cards(:,1));

for i = 1:n
    card_c = deblank(cards(i,:));
    disp(['ineffective tests for ' card_c '...']);
    load([RecsLoc_fp card_c '_fp.mat'],'neft');
    disp(neft);
    disp(['effective tests for ' card_c ' against...']);
    % all effective tests against given control
    eft_t = 0;
    for j = 1:n
        card_s = deblank(cards(j,:));
        % no need for inter-card comparison (except for b6c2/3 against
        % b5...check by hand)
        if i ~= j && card_c(2) == card_s(2)
            load([RecsLoc_fn card_c '/' card_s '_fn.mat'],'eft');
            disp(['  ' card_c '(' card_s ')...']);
            disp(eft);
            
            % ensure we're only keeping track of unique effective tests (if
            % test was effective against more than one subject card, only
            % count it once)
            eft_t = unique([eft_t eft]);
        end
    end
    
    disp(['  overlap between ineffective and effective tests for ' card_c '...'])
    disp(intersect(neft,eft_t));
    disp(' ');
end