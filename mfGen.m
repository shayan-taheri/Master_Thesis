function mfGen(conCard,subCard,refRecsDir,subRecsDir,resultsDir)
%% matched filter routine; gen test
%
% NOTE: we don't want to have to worry about path names here; calling
% routine should take care of all of that (setting/creating data dirs)
%
% conCard: name of control card
% subCard: name of subject card
% refRecsDir: location of conCard's reference record file
% subRecsDir: location of subCard's packed dataset
% resultsDir: where to put results of matched filter tests


%% mfTests variables

% directory delimiter
DD = '/'; 

% length of a record (must be the same as nSp in dilonLoadBinSamples)
lenRec = 175000;
% num of points for all ffts should be the same as lenRec for dilonFFT
nfft = lenRec;
% num of recs to analyse
n = 10000;
% num of recs per iteration
nI = 100;

% constants for aligning records (use with dilonDot)
% load gen refSig (we align w/r/t this)
load([refRecsDir DD conCard],'refRec');
[refSig, refSig_pos] = getRefSig(refRec,'sync');
lenRefSig = length(refSig);
c = repmat(1:nI,lenRefSig,1); %columns
s = repmat([(lenRefSig-1):-1:0]',1,nI); %values to subtract from genPos vector

% open subject dataset for reading (for dilonLoadBinSamples)
cardBin = fopen([subRecsDir DD subCard '.bin']);
negSam = fopen([subRecsDir DD subCard '_neg_s.bin']);


%% define filter outputs
% m0 = struct(...
%    'gen', zeros(1,n), ...
%    'genOrig', zeros(1,n), ...
%    'genNorm', zeros(1,n), ...
%    'genPos', zeros(1,n), ...
%    'vars', char(conCard,subCard,refRecsDir,subRecsDir,resultsDir));

m0 = struct(...
    'gen', zeros(1,n), ...
    'genOrig', zeros(1,n), ...
    'genNorm', zeros(1,n), ...
    'genPos', zeros(1,n), ...
    'recs_norm', zeros(1,n), ... % Shayan's Added !
    'recs_AbsMean', zeros(1,n), ... % Shayan's Added !
    'vars', char(conCard,subCard,refRecsDir,subRecsDir,resultsDir));


%% perform tests
for i = 1:nI:n
    tB_0 = clock;
    % get first and last record to be loaded this iteration
    fRec = i; lRec = i+nI-1; 
    fRec_str = num2str(fRec); lRec_str = num2str(lRec);
    disp([' mf tests for records ' fRec_str ':' lRec_str '...']);
    disp('  loading Records...');
    recs = dilonLoadBinSamples(cardBin,negSam,fRec,lRec);
    disp(sprintf('  ...loading time: %0.2f sec', etime(clock,tB_0)));
    
    
    % generic mf test (find max output within 100 points of expected pos)
    tB_1 = clock;
    disp('  generic test...');
    [m0.genOrig(fRec:lRec), m0.genPos(fRec:lRec)] = ...
        dilonMax(dilonFFTFilt(refSig,recs,nfft),refSig_pos-40,refSig_pos+40);
    disp(sprintf('  ...test time: %0.2f sec', etime(clock,tB_1)));

    
    % extract aligned records
    tB_2 = clock;
    disp(['  aligning records...']);
    p = repmat(m0.genPos(fRec:lRec),lenRefSig,1);
    r = p-s; %rows
    ind = sub2ind([lenRec nI],r,c);
    recs_aligned = recs(ind); %keep these around to use later
    disp(sprintf('  ...aligning time: %0.2f sec', etime(clock,tB_2)));
    
    % perform dotp-based mf test (to ensure that m0.gen == m0.genOrig); get
    % normalised gen mf output, too
    tB_3 = clock;
    disp(['  generic test (dot product)...']);
    [m0.gen(fRec:lRec), m0.genNorm(fRec:lRec)] = ...
        dilonDot(refSig,recs_aligned);
    disp(sprintf('  ...test time: %0.2f sec', etime(clock,tB_3)));
        
    disp(sprintf(' ...iteration time: %0.2f sec', etime(clock,tB_0)));
    
    m0.recs_norm(:,fRec:lRec) = sqrt(sum(recs_aligned.^2,1)); % Shayan's Added !
    
    m0.recs_AbsMean(:,fRec:lRec) = mean(abs(recs_aligned)); % Shayan's Added !
    
end

% save m0 struct
save([resultsDir DD conCard '_' subCard],'m0');

% close files
fclose(cardBin); fclose(negSam);

function [refSig, refSig_pos] = getRefSig(refRec,refSigType)
% load refSig specified by refSigType from refRecsDir

if strcmp('syncAvg',refSigType)
    refSig = wrev(refRec.avg);
    refSig_pos = refRec.sync(2);
    return;
end

% determine start and end points of refSig specified by refSigType
if strcmp('sync',refSigType)
    fSp = refRec.sync(1); lSp = refRec.sync(2);
elseif strcmp('trans1',refSigType)
    fSp = refRec.trans1(1); lSp = refRec.trans1(2);
elseif strcmp('trans1_sync',refSigType)
    fSp = refRec.trans1(1); lSp = refRec.sync(2);
elseif strcmp('trans1_sync_trans2',refSigType)
    fSp = refRec.trans1(1); lSp = refRec.trans2(2);
elseif strcmp('trans1_sync_trans2_mac',refSigType)
    fSp = refRec.trans1(1); lSp = refRec.mac(2);
elseif strcmp('sync_trans2',refSigType)
    fSp = refRec.sync(1); lSp = refRec.trans2(2);
elseif strcmp('sync_trans2_mac',refSigType)
    fSp = refRec.sync(1); lSp = refRec.mac(2);
elseif strcmp('trans2',refSigType)
    fSp = refRec.trans2(1); lSp = refRec.trans2(2);
elseif strcmp('trans2_mac',refSigType)
    fSp = refRec.trans2(1); lSp = refRec.mac(2);
elseif strcmp('mac',refSigType)
    fSp = refRec.mac(1); lSp = refRec.mac(2);
else
    disp('INVALID refSigType! using noise.')
    fSp = 1; lSp = 2000;
end

% build our refSig
refSig = wrev(refRec.rec(fSp:lSp,1));
refSig_pos = lSp;
