function [m0, m0_norm] = dilonDot(refSig,recs)
% dot-product based matched filter; take norm of data, too
% recs: records to perform mf operation on; aligned records assumed
% refSig: refSig used to determine alignment of recs (assumed to be
% reversed and a row vector)

% reverse refSig and transpose (should be row vector)
refSig = wrev(refSig)';

% calc. euclidean norm
refSig_norm = norm(refSig);
recs_norm = sqrt(sum(recs.^2));
% perform dot-product (dotp = A*x)
m0 = refSig*recs;
% normalise output
m0_norm = m0./(refSig_norm*recs_norm);