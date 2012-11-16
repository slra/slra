clear all; d1 = 3; d2 = 3; ell = 2; 
p1 = conv([4 2 1], [5 2]) + [0.05 0.03 0.04 0];
p2 = conv([4 2 1], [5 1]) + [0.04 0.02 0.01 0];
s.m = [d2 - ell + 1; d1 - ell + 1]; 
s.n = d1 + d2 - ell + 1;
z1 = zeros(d2 - ell, 1); z2 = zeros(d1 - ell, 1); 
p = [z1; p1(:); z1; z2; p2(:); z2];
s.w = 1 ./ p; s.w(~isinf(s.w)) = 1;
r = d1 + d2 - 2 * ell + 1;
[ph, info] = slra(p, s, r);
ph1 = ph(2:5), ph2 = ph(8:11), r_ph1 = roots(ph1), r_ph2 = roots(ph2)
d1 = 2; d2 = 2; ell = 1; 
p1 = conv([1   -1], [5   -1]);
p2 = conv([1.1 -1], [5.2 -1]);
s.m = [d2 - ell + 1; d1 - ell + 1]; 
s.n = d1 + d2 - ell + 1;
z1 = zeros(d2 - ell, 1); z2 = zeros(d1 - ell, 1); 
p = [z1; p1(:); z1; z2; p2(:); z2];
s.w = 1 ./ p; s.w(~isinf(s.w)) = 1;
r = d1 + d2 - 2 * ell + 1;
[ph, info] = slra(p, s, r);
ph1 = ph(2:4), ph2 = ph(7:9), r_ph1 = roots(ph1), r_ph2 = roots(ph2)
