more off
N = 2000
c = rand(1, 16) * pi() * 15/16 + pi()/32;
d = (1:N)';
p = d * c;
f0 = sin(d * c) * ones(16,1);
fn = normrnd(0, 0.1, N, 1);


f = f0 +fn;
rk = 32;
Ls = rk+1;
Ks = N-Ls+1;

h0 = hankel(f0(1:Ks),f0(Ks:N));
h = hankel(f(1:Ks),f(Ks:N));
svd(h0)
svd(h)



L= 0.2 *N;
K = N-L+1;

hbig = hankel(f(1:K), f(K:N));
[U,S,V] = svd(hbig);
hrec = U(:, 1:rk) * S(1:rk,1:rk) * V(:, 1:rk)';

frec = zeros(N, 1);
wrec = zeros(N, 1);   

for i = 1:L
	frec(i:(i+K-1),:) = frec(i:(i+K-1),:) + hrec(:,i);
	wrec(i:(i+K-1),:) = wrec(i:(i+K-1),:) + ones(K,1);
endfor
frec = frec ./ wrec;

h1 = hankel(frec(1:Ks,1),frec(Ks:N,1));


x0 = stls(h0(:,1:rk), h0(:,(rk+1):Ls), [2 Ls 1]);
x1 = tls(h1(:,1:rk), h1(:,(rk+1):Ls));


opt.maxiter = 600;
opt.epsgrad = 1e-6;
opt.epsrel = 1e-6;
xs = stls(h(:,1:rk), h(:,(rk+1):Ls), [2 Ls 1], x1,opt)




