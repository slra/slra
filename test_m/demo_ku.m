addpath support

N = 119
rk = 2
Ls = rk+1
Ks = N - Ls+1

f0 = sin((1:N)/12*pi);
fn = normrnd(0, 0.25, 1, N);
f = f0 + fn;
h0 = hankel(f0(1:Ks), f0(Ks:N));
h = hankel(f(1:Ks), f(Ks:N));
[U,S,V0] = svd(h0);
x0 = V0(:,Ls);
x0 = -(x0(1:(Ls-1))/x0(Ls))

a = h(:,1:(Ls-1));
b = h(:,Ls);

[U,S,V] = svd(h);
x = V(:,Ls);
x = -(x(1:(Ls-1))/x(Ls))

s = [2 Ls 1];
pause


h0(1:4,:)
h(1:4,:)

[xh, info, v] = stls(a, b, s, x0);
[ dp, ch, ph ] = corr( xh, h, s);


f_stls = (f'-dp)'
plot([f0' f' f_stls']);

L= floor((N+1)/2)
K = N-L+1

hbig = hankel(f(1:K), f(K:N));
[U,S,V] = svd(hbig);
hrec = U(:, 1:rk) * S(1:rk,1:rk) * V(:, 1:rk)';

size(hbig)



