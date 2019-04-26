% pack records as binary file; trim them in the process

% recsDir: location of unpacked records
% binDir: location to save binary file to

% LOCATION of input/output data
% location of unpacked records
recsDir = '/Volumes/stg/research/2010-SPRING-01';
% location to save binary file to
binDir = '/Volumes/stg/research/2010-SPRING-01';

% num of cards for each model
% b4N = 6; b5N = 10; b6N = 4;
b4N = 6; b5N = 11; b6N = 5;
cardsN = b4N + b5N + b6N; %total num of cards

t0 = clock; %time to run all tests
disp(['Aligning packed recs for ' num2str(cardsN) ' cards...']);

for i = 1:1%cardsN %iterate over cards
    t1 = clock;
    
    % determine card based on counter
    switch logical(true)
        case i <= b4N
            card = ['b4c' num2str(i)];
        case i > b4N && i <= (b4N + b5N)
            card = ['b5c' num2str(i-b4N)];
        case i > (b4N + b5N) && i <= cardsN
            card = ['b6c' num2str(i-b5N-b4N)];
        otherwise
            disp('invalid card, exiting!'); break;
    end
    
    % run conversion routine on recs
    disp(['recs2bin on ' card '...']);
    recs2bin(card,recsDir,binDir);
    disp('...finished');
end