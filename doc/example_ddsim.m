clear all, close all, randn('seed', 0), rand('seed', 0)
Q0 = [1 -1 1]; P0 = [0.81 -1.456 1];  ell = 2; m = 1; q = 2; T = 100; T2 = 40; nl = 0.05;
sys0 = tf(fliplr(Q0), fliplr(P0), 1); 
u0 = rand(T, m); y0 = lsim(sys0, u0); E = rand(T, q); 
w  = [u0 y0] + nl * E / norm(E, 'fro') * norm([u0 y0], 'fro');
u2 = [zeros(ell, 1); ones(T2 - ell, 1)];
y2 = [zeros(ell, 1); NaN * ones(T2 - ell, 1)]; w2 = [u2 y2];
s.m = [ell + 1, ell + 1]; s.n = [T - ell, T2 - ell]; 
s.w = [ones(q * T, 1); inf * vec(~isnan(w2))];  
opt.solver = 'm'; [wh, info] = slra([vec(w); vec(w2)], s, q * ell + 1, opt);
sh = wh(T * q + T2 * m + ell + 1:end); 
s0 = step(sys0, T2 - ell - 1);
plot(sh(2:end), '--b'), hold on, plot(s0(2:end), '-r')
ax = axis; axis([1, T2 - ell - 1, ax(3:4)]), print_fig('slra-ext-f2')
% test_overview_ddsim.m
