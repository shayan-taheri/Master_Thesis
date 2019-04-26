function getDataCont(src,event)
% grab data from previous scan, average and store it, then iniate a new
% scan
global temp time i;

temp(i,:) = mean(event.Data);
time(i,1) = datenum(clock);
i=i+1;

