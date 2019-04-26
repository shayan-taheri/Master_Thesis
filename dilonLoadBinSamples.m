function recs = dilonLoadBinSamples(cardBin,negSam,fRec,lRec)
% load recs from packed binary file; convert to double
% cardBin: fid of packed binary file
% negSam: fid of which samples are negative
% fRec: first rec to load
% lRec: last rec to load

% num of bytes per sample point
nbytes = 1; %sample points stored as uint8
% data type of sample points
dtype = 'uint8';
% num of sample points per rec
nSp = 55000;
% num of recs to return
nRecs = lRec-fRec+1;
% the y-increment used on the scope the data was taken on (need in order to
% convert from int16 to double)
yinc = 0.0169;
% recs to return
recs = zeros(nSp,nRecs);

% determine position of fRec in binary file
offset = nbytes*nSp*(fRec-1);
% move fpos to beginning of fRec
fseek(cardBin,offset,'bof');

% read recs and convert back to original voltage values
recs = yinc*fread(cardBin,[nSp nRecs],dtype);

% restore negative values to rec
offset = nSp/8*(fRec-1); %nSp had better be evenly divisible by 8...fseek only works with bytes
fseek(negSam,offset,'bof');

neg_s = logical(fread(negSam,[nSp,nRecs],'ubit1'));
recs(neg_s) = -recs(neg_s);

% pad recs with random data if dataset is truncated (b4c6)
if size(recs,2) ~= nRecs
    disp('WARNING: TRUNCATED DATASET! Appending random data.')
    recs = [recs, rand(nSp,nRecs-size(recs,2))];
end