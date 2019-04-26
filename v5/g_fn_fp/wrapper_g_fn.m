% call wrapper function to find fn for all cards
% NOTE: fp must be known for cards and 'matlabpool(2)' should be run at the
% command window before proceeding

% call wrapper function to find fp for all cards
% P = eval('pwd');
% make base files available
% addpath([P '/../base']);

addpath('/home/shayan/pli/v5/base');

% LOCATION of input/output data
% location of filter output data
subDir = '/home/shayan/pli/data/SP1/tests_results/1';
% dir for results of fp tests (reading)
fpDir = '/home/shayan/pli/data/SP1/fn_fp/fp_test';
% dir for results of fn tests (writing)
fnDir = '/home/shayan/pli/data/SP1/fn_fp/fn_test';
% directory delimiter
DD = '/'; 

% num of cards for each model; assume same number of cards for each dataset
b4N = 6; b5N = 10; b6N = 10;
cardsN = b4N + b5N + b6N; %total num of cards


for i = 1:cardsN %iterate over subject card
    
    % determine control card based on counter
    switch logical(true)
        case i <= b4N
            conCard = ['b4c' num2str(i)];
        case i > b4N && i <= (b4N + b5N)
            conCard = ['b5c' num2str(i-b4N)];
        case i > (b4N + b5N) && i <= cardsN
            conCard = ['b6c' num2str(i-b5N-b4N)];
    end
    
    % load fp data for card
    load([fpDir DD conCard '_' conCard],'th');

    for j = 1:cardsN %iterate over control card
        % determine control card based on counter
        switch logical(true)
            case j <= b4N
                subCard = ['b4c' num2str(j)];
            case j > b4N && j <= (b4N + b5N)
                subCard = ['b5c' num2str(j-b4N)];
            case j > (b4N + b5N) && j <= cardsN
                subCard = ['b6c' num2str(j-b5N-b4N)];
        end
	% comment out if running inter-dataset comparisons
%        if i ~= j
	    disp([conCard ' vs. ' subCard '...']);
            g_fn(conCard, subCard, subDir, fnDir, th);
%        end
    end
%     % no need for inter-card comparison (except for b6c2/3 against b5...check by hand)
%     for j = 1:n
%         card_s = deblank(cards(j,:));
%         if i ~= j && card_c(2) == card_s(2)      
%             disp('Finding FN for ...');
%             wrapper_fn_func(card_c,card_s,th,fu_dfpi.th,[]);
%         end
%     end
end
