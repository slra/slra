% Test of the program of Nicola

clear all
n  = 2;
np = 10;
M  = round(linspace(50,500,np));

% true data
m = M(end);
z = linspace(0.9,1,n);
v = zeros(m,n);
for i = 1:m
  v(i,:) = z.^(i-1);
end
h = v * diag(ones(n,1)) * v(1:n+1,:)';

a0 = h(:,1:n);
b0 = h(:,n+1);
c0 = fliplr([a0 b0]);
a0 = c0(:,1:n);
b0 = c0(:,n+1);
x0  = a0\b0;
nx0 = norm(x0);
p0  = [c0(1,end:-1:2)'; c0(:,1)];

% perturbed data
pt = .05 * rand(m+n,1);
p  = p0 + pt;
c  = toeplitz(p(n+1:end),p(n+1:-1:1));
a  = c(:,1:n);
b  = c(:,n+1:end);

t = zeros(1,np);
f = zeros(1,np);
e = zeros(1,np);
for i = 1:np
  i = i
  % estimation
  xtls = tls(a,b);
  flops(0), 
  tic, 
  x = stln_s(a(1:M(i),:),b(1:M(i)),xtls);
  t(i) = toc;
  f(i) = flops;
  e(i) = norm(x0 - x) / nx0;  
end

[M' f']
plot(M,f)


