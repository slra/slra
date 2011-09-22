addpath ../matlab_version

m = 50000;
n = 100;
d = 10;

a = rand(m,n);
b = rand(m,d);

tic, x1 = tls(a,b); t1 = toc
tic, x2 = efftls(a,b); t2 = toc
tic, x3 = mb02md([a b],n); t3 = toc

norm(x1-x2,'fro')
