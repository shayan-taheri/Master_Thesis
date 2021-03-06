% display ineffecitve tests, all cards

% DATA INFORMATION
RecsLoc = '/Users/gerdes/dilon/';
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
    load([RecsLoc_fp card_c '_fp.mat'],'neft');
    disp(['ineffective tests for ' card_c '...']);
    disp(neft);
end