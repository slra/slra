clear all; randn('seed', 0); load aut_sys_traj; T = length(y0); ell = 8; 
plot(y0(1:T), 'r-'), ax = axis; axis([0 T ax(3:4)]), print_fig('slra-f1')
size(null(blkhank(y0(1:T), ell + 1)')', 1) 
y2R = @(y, ell) null(blkhank(y, ell + 1)')'; R0 = y2R(y0, ell);
y = y0 + 0.3 * std(y0) * randn(length(y0), 1); 
s.m = ell + 1; s.n = T - ell; r = ell; 
y2R_appr = @(y, ell) fliplr(poly(eig(h2ss(y, ell))));, opt.Rini = y2R_appr(y(1:T), ell);
[yh, info] = slra(y(1:T), s, r, opt); 
plot(yh, 'b--'), hold on, plot(y(1:T), 'k:'), plot(y0(1:T), 'r-')
axis([1 T ax(3:4)]), print_fig('slra-f2')
