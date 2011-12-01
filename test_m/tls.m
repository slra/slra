% TLS --- TLS solution of Ax = B

function [xh,ch] = tls(a,b)

n  = size(a,2);
d  = size(b,2);

[u,s,v]  = svd([a,b],'econ');
xh = -v(1:n,n+1:end) / v(n+1:end,n+1:end);

if nargout == 2
  ch = u(:,1:n) * s(1:n,1:n) * v(:,1:n)';
end
