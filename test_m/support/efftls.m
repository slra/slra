% TLS --- TLS solution of Ax = B; efficinet version

function xh = tls(a,b)

n  = size(a,2);
d  = size(b,2);

r = triu(qr([a b]));
a = r(1:n+d,1:n);
b = r(1:n+d,n+1:end);

[u,s,v]  = svd([a,b],0);
xh = -v(1:n,n+1:end) / v(n+1:end,n+1:end);
