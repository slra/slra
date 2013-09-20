clear all, n1 = 2; n2 = 2; n3 = 2; ell = 1; 
p1 = conv([1   -1], [5   -1]);
p2 = conv([1.1 -1], [5.2 -1]);
p3 = conv([1.2 -1], [5.4 -1]);
s.m = [n1 - ell + 1; n1 + n2 + n3 - 2 * ell + 2]; 
s.n = 2 * n1 + 2 + n3 - 2 * ell + 2;
s.phi = eye(sum(s.m));
s.phi((s.m(1) + n2 - ell + 2):(s.m(1) + n2 - ell + 1 + n1), :) = [];
z1 = inf * ones(n1 - ell, 1); z2 = inf * ones(s.m(2) - 1, 1);
s.w = [z1; ones(n2 + 1, 1); z1; ones(n3 + 1, 1); 
       z1; z2; ones(n1 + 1, 1); z2];
r = s.m(1) + n2  + n3 - 3 * ell + 2;
p = zeros(size(s.w)); p(s.w == 1) = [p1 p2 p3]';
opt.solver = 'r'; [ph, info] = slra(p, s, r, opt);
ph123 = ph(s.w == 1); 
r_ph1 = roots(ph123(1:n1 + 1))
r_ph2 = roots(ph123(n1 + 2:n1 + n2 + 2))
r_ph3 = roots(ph123(n1 + n2 + 3:end))
