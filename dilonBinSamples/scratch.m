cd /local/rgerdes
fid=fopen('b4c1.bin','w')

open b4c1/sample00001.mat
rec = ans.rec'; %why did we store record as column vector? 
                %Look to analysis code to find correct form (row vs column)
                %fftfilt applies filter to columns, so we store as column
                %vector
rec = rec(1:55000);
yinc=0.02
rec_int16 = int16(rec/yinc);

% want to store samples using as few bits as possible; can use uint8 if we
% get rid of negative values, so we remember which samples were negative
n_samples=(rec_int16<0);

rec_uint8=uint8(abs(rec_int16));

fwrite(fid,rec_uint8,'uint8');
fwrite(fid,n_samples,'ubit1');
fclose(fid);

% read things back
fid=fopen('b4c1.bin','r');
rec8 = fread(fid,[length(rec),1],'uint8'); %matlab appears to convert to double by default
n_s = logical(fread(fid,[length(rec),1],'ubit1'));
r = rec8;
r(n_s) = -r(n_s);
r = yinc*r;
fclose(fid);