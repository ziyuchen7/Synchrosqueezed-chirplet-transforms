% imageSQ.m
%
% Display time-frequency result of the Synchrosqueezing transform and others
% function imageSQ(t, ytic, M) ;
%
function imageSQ(t, ytic, M, d) ;

fz = 22;

S = size(M);
Q = M(:);
q = quantile(Q, d);
M(find(M>q)) = q;

imagesc(t, ytic, M)
axis xy ;
set(gca, 'fontsize', fz);
