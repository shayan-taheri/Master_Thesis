% get fnr for all all cards

% LOCATION of input/output data
% location of filter output data
% dir for results of fn tests
fnDir = '/Users/rgerdes/Documents/dilon/research/eval_2_2/v5/results/2010-SPRING-01/1';
% directory delimiter
DD = '/'; 

% num of cards for each model; assume same number of cards for each dataset
b4N = 6; b5N = 11; b6N = 10;
cardsN = b4N + b5N + b6N; %total num of cards

% tests to calculate fnr for
% tsts = [1    60    61    62    63   329   330   331   332]; %b4
% tsts = [1   255   256   257   258   259   288]; %b5
tsts = [1    60    61    62    63   329]; %b6

% model to use of intra-model comparison
m = 'b6';

fid = fopen('fnr.txt','w');
for i = 1:cardsN %iterate over control card
    
    % determine control card based on counter
    switch logical(true)
        case i <= b4N
            conCard = ['b4c' num2str(i)];
        case i > b4N && i <= (b4N + b5N)
            conCard = ['b5c' num2str(i-b4N)];
        case i > (b4N + b5N) && i <= cardsN
            conCard = ['b6c' num2str(i-b5N-b4N)];
    end
    
    for j = 1:cardsN %iterate over subject card
        % determine subject card based on counter
        switch logical(true)
            case j <= b4N
                subCard = ['b4c' num2str(j)];
            case j > b4N && j <= (b4N + b5N)
                subCard = ['b5c' num2str(j-b4N)];
            case j > (b4N + b5N) && j <= cardsN
                subCard = ['b6c' num2str(j-b5N-b4N)];
        end
        
        if i ~= j && strcmp(subCard(1:2),m) && strcmp(conCard(1:2),m)
            load([fnDir DD conCard '_' subCard],'acpt_cnt','ns');
            % calc fnr
            fnr = num2str(100*mean(acpt_cnt(:,tsts))/ns);
            fnr_str = [conCard ' vs. ' subCard ' fnr: ' fnr '%'];
            disp(fnr_str);
            % write fnr to file
            fprintf(fid,'%s\n',fnr_str);
        end
    end
    disp(' ');
    fprintf(fid,'\n');
end
fclose(fid);