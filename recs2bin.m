% function rec_ret = recs2bin(card,recsDir,binDir)
% convert dilon recs to a packed binary file (convert to uint8 first)
% card: name of card
% recsDir: location of unpacked records
% binDir: location to save binary file to

card = 'b4c5';
recsDir = '/home/shayan/pli/data/10mb/2013-SPRING-01/b4/b4c5';
binDir = '/home/shayan/pli/data/SP1/pack_recs';

% recs to convert and pack
fRec = 1; lRec = 10000;
% sample points to include in packed record
fSp = 1; lSp = 55000; %NOTE: this needs to be evenly divisible by 8, otherwise our read routine won't work
% the y-increment used on the scope the data was taken on (need in order to
% convert from double to integer)
 yinc = 0.02;
% how often to update status, every n-records
n = 1000;
tRec = lRec-fRec+1; tRecStr = num2str(tRec);

cardBin = fopen([binDir '/' card '.bin'],'w+');
negSam = fopen([binDir '/' card '_neg_s.bin'],'w+');

for i = fRec:lRec
    if mod(i,n) == 0
        disp([card ' ' num2str(i/tRec*100) '% complete (' num2str(i) '/' tRecStr ' recs converted)...']);
    end
    

    % ************************* Shayan's Modification *************************
    
    %{
    % load each record
    iStr = num2str(i);
    recStr = [recsDir '/sample' ...
        strrep(num2str(zeros(1,5-length(iStr))), ' ', '') iStr '.mat'];
    load(recStr,'rec');
    %}
    
    % load([recsDir,'/sample','00001','.mat'],'rec');
    
    if (i >= 1) && (i <= 9)
        
        load([recsDir,'/sample','0000',num2str(i),'.mat'],'rec');
        
    elseif (i >= 10) && (i <= 99)
        
        load([recsDir,'/sample','000',num2str(i),'.mat'],'rec');
        
    elseif (i >= 100) && (i <= 999)
        
        load([recsDir,'/sample','00',num2str(i),'.mat'],'rec');
        
    elseif (i >= 1000) && (i <= 9999)
        
        load([recsDir,'/sample','0',num2str(i),'.mat'],'rec');
    
    else
        
        load([recsDir,'/sample',num2str(i),'.mat'],'rec');
        
    end

    % ************************* Shayan's Modification *************************


    % only work with samples we're going to write to bin file
    rec = rec(fSp:lSp)'; %transpose rec so that it's a row vector (for fftfilt)

    % want to store samples using as few bits as possible; can use uint8 if we
    % get rid of negative values, so we remember which samples were negative
    neg_s=(rec < 0);
    
    % removoe neg vals, make samples integer, and convert to uint8
    rec_uint8=uint8(abs(rec/yinc));
    
    % write rec as uint8 binary data
    fwrite(cardBin,rec_uint8,'uint8');

    % write location of neg samples
    fwrite(negSam,neg_s,'ubit1');
end

fclose(cardBin);
fclose(negSam);
