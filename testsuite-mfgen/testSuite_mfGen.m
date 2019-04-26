% testSuite_mfGen.m: based upon testSuite.m v5; only perform gen mf test

% LOCATION of input/output data: individual datasets occupy a single row
% location of packed binary data records
recsDir = ['/home/shayan/pli/data/SP1/pack_recs'];
% reference recs location; NOTE: must be hand selected for each test
refRecsDir = ['/home/shayan/pli/data/SP1/ref_recs'];
% dir for results of tests
resultsDir = ['/home/shayan/pli/data/SP1/test_res'];
% directory delimiter
DD = '/'; 

% num of cards for each model; assume same number of cards for each dataset
b4N = 6; b5N = 10; b6N = 10;
cardsN = b4N + b5N + b6N; %total num of cards

% num of comparisons (control vs. subject)
nRecsDir = size(recsDir,1);
nComp = (nRecsDir-1)*cardsN^2;
nCompStr = num2str(nComp);
nCompFin = 0; % num of comparisons made in previous runs

t0 = clock; % time to run all tests
disp(['mf testsuite (asymmetric) for ' num2str(cardsN) ' cards and ' ...
    num2str(size(recsDir,1)) ' datasets totals ' nCompStr ' comparisons...']);

% routine overview: iterate over datasets and select control card
% select subject card (parallel execution) and iterate over datasets

% iterate over datasets (i.e. use refRecs from one dataset against self and
% other dataset).  if we only want to run first dataset against the rest,
% set i = 1:1
for i = 1:size(recsDir,1) %iterate over datasets
    % select control card
    for j = 1:cardsN %iterate over cards
        
        % determine control card based on counter
        switch logical(true)
            case j <= b4N
                conCard = ['b4c' num2str(j)];
            case j > b4N && j <= (b4N + b5N)
                conCard = ['b5c' num2str(j-b4N)];
            case j > (b4N + b5N) && j <= cardsN
                conCard = ['b6c' num2str(j-b5N-b4N)];
        end
        
        % select subject card and run tests against all data sets.
        % NOTE ON PARALLEL IMPLEMENTATION: we iterate over subject cards
        % using parfor so that no two instances of matlab are attempting to
        % access a dataset for the same card in the hopes of avoiding  I/O 
        % contention
        for k = 1:cardsN % iterate over cards
%         parfor k = 1:cardsN %iterate over cards
            switch logical(true)
                case k <= b4N
                    subCard = ['b4c' num2str(k)];
                case k > b4N && k <= (b4N + b5N)
                    subCard = ['b5c' num2str(k-b4N)];
                case k > (b4N + b5N) && k <= cardsN
                    subCard = ['b6c' num2str(k-b5N-b4N)];
            end
            
            % finally iterate over datasets for subject
            for l = 1:size(recsDir,1) % iterate over datasets (l = 2:2)
                % in case things go wrong and we need to restart, remember
                % i, j, k,l
                disp(['(debug) loop positions (i j k l): ' num2str(i) ' ' num2str(j) ...
                    ' ' num2str(k) ' ' num2str(l)]);
                % which cards, and from which datasets, are control and
                % subject
                disp(['(status) control: ' conCard ' of ' recsDir(i,:)]);
                disp(['(status) subject: ' subCard ' of ' recsDir(l,:)]);
                
%                mfGen(conCard,subCard, ...
%                    refRecsDir(i,:), recsDir(l,:), ...
%                    [resultsDir(i,:) DD num2str(l+2)])
                
                 mfGen(conCard,subCard, ...
                     refRecsDir(i,:), recsDir(l,:), ...
                     [resultsDir(i,:) DD num2str(l)])

                % where are we in this mess...has to be after parfor
                % because parfor executes out of order; uncomment for
                % single loop
                nCompMade = (i-1)*cardsN^2*(nRecsDir-1) + ...
                    (j-1)*cardsN*(nRecsDir-1) + ...
                    (k-1)*(nRecsDir-1) + 1; % comparisons made thus far
                disp(['(status) ' num2str(nCompMade) '/' nCompStr ...
                    ' (' num2str(nCompMade/nComp*100) '%)'])
                t1 = etime(clock,t0);
                disp(sprintf('(status) running time: %0.2f sec', t1));
                % estimate time remaining (avg. time per test * remaining
                % tests)
                est = t1/(nCompMade-nCompFin)*(nComp-nCompMade);
                disp(sprintf('(status) estimated time to completion: %0.2f sec', est));
                disp(' '); disp(' ');
            end
        end
        
%         % where are we in this mess...has to be after parfor
%         % because parfor executes out of order
%         nCompMade = (i-1)*cardsN^2*nRecsDir + ...
%             (j-1)*cardsN*nRecsDir + ...
%             cardsN*nRecsDir;%comparisons made thus far
%         disp(['(status) ' num2str(nCompMade) '/' nCompStr ...
%             ' (' num2str(nCompMade/nComp*100) '%) comparisons made'])
%         t1 = etime(clock,t0);
%         disp(sprintf('(status) running time: %0.2f sec', t1));
%         % estimate time remaining (avg. time per test * remaining
%         % tests)
%         est = t1/nCompMade*(nComp-nCompMade);
%         disp(sprintf('(status) estimated time to completion: %0.2f sec', est));
    end
end
