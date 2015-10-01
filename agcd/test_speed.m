ntest = 8;
times = zeros(ntest,3);
iters = zeros(ntest,3);
res = zeros(ntest,3);

h0 = [-1; 1; 0; 10; 1];
d = length(h0) - 1;

times(:,1) = 1:ntest;
iters(:,1) = 1:ntest;
res(:,1) = 1:ntest;
for k=1:ntest
  k
  n1 = 25*k; n2 = 15*k; n3 = 10*k;
  g1 = conv(conv([-1; zeros(n1 - 1, 1); 1], ...
      [-2; zeros(n2 - 1, 1) ; 1]), [-3; zeros(n3 - 1, 1) ; 1]);
  g2 = conv(conv([1; zeros(n1 - 1, 1); 1], ...
      [5; zeros(n2 - 1,1) ; 1]), [2; zeros(n3 - 1, 1) ; 1]);
  p = conv(g1, h0); q = conv(g2, h0);
  sigma = 1e-2;
  p = p + sigma * randn(size(p));
  q = q + sigma * randn(size(q));
  
  
  gini = g_ini({p,q}, d);
  h = lsdivmult({p,q}, d, gini)';
  g1 = gini(1:length(p)-d)';
  g2 = gini(length(p)-d + (1:length(q)-d))';
  hh = h; g1h = g1; g2h = g2;

  opt.hini = hh';
  opt.epsgrad = 1e-20;
  opt.maxiter = 1;

  [ph, info] = gcd_nls({p,q}, [], d, opt);
  res(k,2) =  norm([p;q] -cell2mat(ph),2);
  iters(k,2) = info.iter;
  
  [hh, g1h, g2h, resn, iter] = c_f_newton_iter_mod(p',q', ...
                                   hh,g1h,g2h, 0, opt.maxiter+1, 0);
  iters(k,3) = iter-1;
  res(k,3) =  norm([p;q] -[conv(hh,g1h)';conv(hh,g2h)'],2);

  hh = h; g1h = g1; g2h = g2;
  times(k,2) = meas_func(@gcd_nls, {{p,q}, [], d, opt});
  times(k,3) = meas_func(@c_f_newton_iter_mod, {p',q', ...
                               hh,g1h,g2h, 0, opt.maxiter+1, 0});
end

save('res_speed.txt', 'res', '-ascii');
save('iter_speed.txt', 'iters', '-ascii');
save('times_speed.txt', 'times', '-ascii');