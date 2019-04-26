% display effecitve tests, all cards

% DATA INFORMATION
RecsLoc = '/Users/gerdes/dilon/';
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

for i = 1:n
    card_c = deblank(cards(i,:));
    disp(['effective tests for ' card_c ' against...']);
    % no need for inter-card comparison (except for b6c2/3 against b5...check by hand)
    for j = 1:n
        card_s = deblank(cards(j,:));
        if i ~= j && card_c(2) == card_s(2)
            load([RecsLoc_fn card_c '/' card_s '_fn.mat'],'eft');
            disp(['  ' card_c '(' card_s ')...']);
            disp(eft);
        end
    end
end