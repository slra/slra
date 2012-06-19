clear all; randn('seed', 0); load aut_sys_traj; T = length(y0); ell = 8; 
plot(y0(1:T), 'r-'), ax = axis; axis([0 T ax(3:4)]), print_fig('slra-f1')
size(null(blkhank(y0(1:T), ell + 1)')', 1) 
y2R = @(y, ell) null(blkhank(y, ell + 1)')'; R0 = y2R(y0, ell);
y = y0 + 0.25 * std(y0) * randn(length(y0), 1); 
s.m = ell + 1; s.n = T - ell; r = ell; 
y2R_appr = @(y, ell) fliplr(poly(eig(h2ss(y, ell))));, opt.Rini = y2R_appr(y(1:T), ell);
[yh, info] = slra(y(1:T), s, r, opt); 
plot(yh, 'b--'), hold on, plot(y(1:T), 'k:'), plot(y0(1:T), 'r-')
axis([1 T ax(3:4)]), print_fig('slra-f2')
yc = y(1:T) - mean(y(1:T)); opt_c.Rini = fliplr(poly(eig(h2ss(yc, ell))));
[yhc, info_c] = slra(yc, s, r, opt_c); info_c.fmin
ell = ell + 1; s.m = ell + 1; s.n = T - ell; r = ell; 
zf = 1; opt_f.psi = sylvester(poly(zf), 1, ell - length(zf) + 1); 
opt_f.Rini = opt_c.Rini * opt_f.psi;
[yh_f, info_f] = slra(y(1:T), s, ell, opt_f); info_f.fmin
ell = 8; s.m = ell + 1; s.n = T - ell; r = ell; 
nt = round((ell + 1) / 2); I = eye(nt); 
opt_p.psi = [I flipud(I(:, (1 + mod(ell + 1, 2)):end))];
opt_p.Rini = fliplr(poly(eig(approx(y(1:T), ell))));
[yh_p, info_p] = slra(y(1:T), s, ell, opt_p);
abs(roots(info_p.Rh))
ell = 8; T = [20 30 40]; N = length(T); ind = [1 50 100];
s.m = ell + 1; s.n = T - ell; r = ell;   
p = []; for i = 1:N, p = [p; y(ind(i):ind(i) + T(i) - 1)]; end
[ph_m, info_m] = slra(p, s, r, opt); 
angle = @(a, b) min(abs(acosd(a(:)' * b(:) / norm(a) / norm(b))), ...
                    abs(180 - acosd(a(:)' * b(:) / norm(a) / norm(b))));, angle(info_m.Rh, R0) 
s = rmfield(s, 'n');
for i = 1:N
  p = y(ind(i):ind(i) + T(i) - 1);
  [ph, info] = slra(p, s, r, opt);
  angle(info.Rh, R0) 
end
T = 100; Im = (ell + 1):(ell + 1):T;
p = y(1:100); p0 = y0(1:100); p(Im) = 0; 
s.w = ones(T, 1); s.w(Im) = 1e-5;
[ph_m, info_m] = slra(p, s, r, opt);
p(Im) = NaN;
figure(3), plot(ph_m, 'b--'), hold on, 
plot(p, 'k:'), plot(p0(1:T), 'r-')
ax = axis; 
plot(Im, ax(3) * ones(size(Im)), ...
     'Xk', 'markersize', 15)
axis(ax), print_fig('slra-f3')
