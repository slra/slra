clear all; randn('seed', 0); load aut_sys_traj; opt.solver = 'r';
size(null(blkhank(y0, ell + 1)')', 1) 
y2R = @(y, ell) null(blkhank(y, ell + 1)')'; R0 = y2R(y0, ell);
nl = 0.4; y = y0 + nl * std(y0) * randn(T, 1); 
rank(blkhank(y, ell + 1)) 
y2R_appr = @(y, ell) fliplr(poly(eig(h2ss(y, ell))));, opt.Rini = y2R_appr(y, ell);
s.m = ell + 1; s.n = T - ell; r = ell; 
[yh, info] = slra(y, s, r, opt);
plot(yh, 'b--'), hold on, plot(y, 'k:'), plot(y0, 'r-')
ax = axis; axis([1 T ax(3:4)]), print_fig('slra-f2')
R = info.Rh; yini = yh(end - ell + 1:end); X = - R(1:ell) / R(end); yvh = yini;
                                           for t = 1:Tv
                                             yvh = [yvh; X * yvh(end - ell + 1:end)]; 
                                           end  
                                           yvh = yvh(ell + 1:end);
figure, plot(yvh, 'b--'), hold on, plot(yv(3:end), 'r-')
ax = axis; axis([1 T ax(3:4)]), print_fig('slra-f2v')
zf = [1 1]; opt.psi = mult_mat(poly(zf), 1, ell - length(zf) + 1); 
[yh_f, info_f] = slra(y, s, ell, opt);
R = info_f.Rh; yini = yh_f(end - ell + 1:end); X = - R(1:ell) / R(end); yvh = yini;
                                               for t = 1:Tv
                                                 yvh = [yvh; X * yvh(end - ell + 1:end)]; 
                                               end  
                                               yvh = yvh(ell + 1:end);, plot(yvh, 'k-.')
%[yvh_f, infov_f] = misfit(yv, s, info_f.Rh); 
ax = axis; axis([1 T ax(3:4)]), print_fig('slra-f2v')
s.m = ell + 1; s.n = T - ell; r = ell; 
Im = (ell + 1):(ell + 1):T; opt.psi = [];
y(Im) = NaN; s.w = ones(T, 1); s.w(Im) = 1e-5;
[yh_m, info_m] = slra(y, s, r, opt);
figure(3), plot(yh_m, 'b--'), hold on, 
plot(y, 'k:'), plot(y0, 'r-')
ax = axis; 
plot(Im, ax(3) * ones(size(Im)), ...
     'Xk', 'markersize', 15)
axis(ax), print_fig('slra-f3')
