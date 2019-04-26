function [m,p] = dilonMax(y,fSp,lSp)
% find max of filter output within specified range; return value and
% position (getting odd alignment on b4)
% y: filter output (all shifts) for recs; sample point-by-record
% fSp: sample point to start looking for max
% lSp: last sample point to look for max
% m: maximum filter output between fSp and lSp
% p: position of m (in terms of overall rec length)

% find max value between fSp:lSP
[m,p] = max(y(fSp:lSp,:));
% know position of max val with offset of fSp, put it in terms of overall
% record
p = p + fSp - 1;