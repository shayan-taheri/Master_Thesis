% create a refRec struct for a card: prompt user for location of certain
% parts of signal and then create ref rec struct
card = 'b6c10';
% input: packed binary records for card
recsDir = '/local/rgerdes/2010-SPRING-02';
% output: locatation to write ref rec
refRecsDir = '/local/rgerdes/2010-SPRING-02/refRecs';
% num of recs to load (use for averaging)
nRecs = 25;
% which rec to use as base of refRec
refRecPos = 1;

% load recs
cardBin = fopen([recsDir '/' card '.bin']);
negSam = fopen([recsDir '/' card '_neg_s.bin']);
recs = dilonLoadBinSamples(cardBin,negSam,refRecPos,refRecPos+nRecs-1);

% create halfscreen plot of refRec
fullscreen = get(0,'ScreenSize');
figure('Position',[fullscreen(3)/2 -50 fullscreen(3)/2 fullscreen(4)])
% figure('Position',[0 -50 fullscreen(3)/2 fullscreen(4)])
plot(recs(:,1));


% get locations of portions of interest
def = input('Use defaults ([y]/n)? ','s');

if (strcmp(def,'') || strcmp(def,'y'))
    % defaults
%     % b4
%     trans1_b = 3000;
%     sync_b = 3275;
%     trans2_b = 18285;
%     mac_b = 18755;
%     mac_e = 53785;
%     % b5
%     trans1_b = 2650;
%     sync_b = 2925;
%     trans2_b = 17915;
%     mac_b = 18400;
%     mac_e = 53415;
    % b6
    trans1_b = 3000;
    sync_b = 3300;
    trans2_b = 18305;
    mac_b = 18800;
    mac_e = 53810;
else
    disp('input first sample point of...');
    trans1_b = input(' transient: ');
    sync_b = input(' sync signal: ');
    trans2_b = input(' transient to mac: ');
    mac_b = input(' mac: ');
    mac_e = input(' mac (end): ');
end

% build refRec struct
refRec = struct(...
    'rec', recs(:,1), ...
    'trans1', [trans1_b; sync_b-1], ...
    'sync', [sync_b; trans2_b-1], ...
    'trans2', [trans2_b; mac_b-1], ...
    'mac', [mac_b; mac_e], ...
    'avg', zeros(length(sync_b:trans2_b),1));

% build refRec_avg (refRec built of nRecs averaged)
% constants for aligning records
refSig = wrev(recs(sync_b:trans2_b,1));
lenRS = length(refSig);
lenRec = size(recs,1);
c = repmat(1:nRecs,lenRS,1); %columns
s = repmat([(lenRS-1):-1:0]',1,nRecs); %values to subtract from genPos vector

% get position of maximum alignment
[m0, m0_pos] = ...
    dilonMax(fftfilt(refSig,recs),trans2_b-49,trans2_b+50);
% extracted aligned portion of recs
p = repmat(m0_pos,lenRS,1);
r = p-s; %rows
ind = sub2ind([lenRec nRecs],r,c);
recs_aligned = recs(ind);
% average recs
refRec.avg = mean(recs_aligned,2);

% save refRec struct and clean up
save([refRecsDir '/' card], 'refRec');
fclose(cardBin); fclose(negSam);
close(gcf);
% clear;