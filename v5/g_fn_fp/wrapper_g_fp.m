% call wrapper function to find fp for all cards
% P = eval('pwd');
% make base files available
% addpath([P '/../base']);

addpath('/home/shayan/pli/v5/base');

% LOCATION of input/output data
% location of filter output data
conDir = '/home/shayan/pli/data/SP1/tests_results/1';
% dir for results of tests
fpDir = '/home/shayan/pli/data/SP1/fn_fp/fp_test';
% directory delimiter
DD = '/'; 

% num of cards for each model; assume same number of cards for each dataset
b4N = 6; b5N = 10; b6N = 10;
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
    disp([conCard '...']);
    g_fp(conCard,conDir,fpDir);
end
