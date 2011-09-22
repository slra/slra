% GTLS --- GTLS solution of Ax = B

function xh = gtls(a,b,V)

[m,n]  = size(a);

cholV    = chol(V);
invcholV = inv(cholV);
ab  = [a b] * invcholV;

[u,s,v]  = svd(ab,0);
v = v(:,n+1:end);
xh = invcholV * v;
xh = -xh(1:n,:) / xh(n+1:end,:);





