function ind = badRec(m0)

% % expect all output to fall within two standard deviations
% ind = vec2ind(m0 < (u_m0 - 2*s_m0))
% u_m0 = mean(m0);
% s_m0 = std(m0);

% if m0 is closer to zero than mean it's bad
d2m = abs(m0-mean(m0));
ind = find(d2m > m0);