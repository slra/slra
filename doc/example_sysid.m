clear all, randn('seed', 0); opt.solver = 'c';
Q0 = [1 -1 1]; P0 = [0.81 -1.456 1];  ell = 2; m = 1; q = 2; T = 100; T2 = 40; nl = 0.05;
sys0 = tf(fliplr(Q0), fliplr(P0), 1); 
u0 = rand(T, m); y0 = lsim(sys0, u0); E = rand(T, q); 
w  = [u0 y0] + nl * E / norm(E, 'fro') * norm([u0 y0], 'fro');
s.m = (ell + 1) * ones(q, 1); s.n = T - ell; r = ell * q + m;
[wh1, info] = slra(vec(w), s, q * ell + 1);
Qh = fliplr(info.Rh(1:ell + 1)); Ph = -fliplr(info.Rh(ell + 2:end));
sysh = tf(Qh, Ph, 1);
s0 = step(sys0, T2 - ell); sh = step(sysh, T2 - ell); 
plot(sh(2:end), '--b'), hold on, plot(s0(2:end), '-r')
legend('sh', 'sb  .', 'Location', 'SouthEast')
ax = axis; axis([1, T2 - 1, ax(3:4)]), print_fig('f-sysid')
