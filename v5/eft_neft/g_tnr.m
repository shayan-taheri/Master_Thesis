% get tnr for all all cards

% LOCATION of input/output data
% dir for results of fp tests
fpDir = '/Users/rgerdes/Documents/dilon/research/eval_2_2/v5/results/2010-SPRING-01/1';
% directory delimiter
DD = '/'; 

% num of cards for each model; assume same number of cards for each dataset
b4N = 6; b5N = 11; b6N = 10;
cardsN = b4N + b5N + b6N; %total num of cards



for i = 1:cardsN %iterate over cards
    
    % determine control card based on counter
    switch logical(true)
        case i <= b4N
            conCard = ['b4c' num2str(i)];
        case i > b4N && i <= (b4N + b5N)
            conCard = ['b5c' num2str(i-b4N)];
        case i > (b4N + b5N) && i <= cardsN
            conCard = ['b6c' num2str(i-b5N-b4N)];
    end
    
    load([fpDir DD conCard '_' conCard],'a_t','nc','n','nb');
    disp([conCard ' tnr: ' num2str(a_t(1)/(nc-n))]);
end