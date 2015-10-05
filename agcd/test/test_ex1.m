clear all;
addpath 'uvgcd';

p = 1;
q = 1;

for j=1:10
  x_j = (-1)^j * j/2;
  p = conv(p,[-x_j; 1]);
  q = conv(q,[-x_j+ 10^(-j); 1]);
end

sccoef = norm(p,2);
p = p./sccoef;
q = q./sccoef;

res = zeros(10,7);
iters = zeros(10,4);
res(:,1) = (1:10)';
iters(:,1) = (1:10)';
res(:,3) = [5.17e-1;6.95e-4;1.97e-5;2.89e-6;5.28e-5;...
            2.15e-3;8.34e-2;2.04e0;4.70e1;7.73e2] / sccoef;
            
opt.disp = 'iter';
opt.epsgrad = 1e-20;


for d=1:5
  d
  [ph, info] = gcd_nls({p,q}, [], d, opt);
  res(d,5) =  norm([p;q] -cell2mat(ph),2);
  iters(d,2) = info.iter;
end

for d=7:10
  [ph, info] = gcd_cofe({p,q}, [], d, opt);
  res(d,5) =  norm([p;q] -cell2mat(ph),2);
  iters(d,2) = info.iter;
end

for d=1:10
  [ph, info] = gcd_syl({p,q}, [], d, opt);
  res(d,6) =  norm([p;q] -cell2mat(ph),2);
  iters(d,3) = info.iter;
end


for d=1:10  
  gini = g_ini({p,q}, d);
  hh = lsdivmult({p,q}, d, gini)';
  g1h = gini(1:length(p)-d)';
  g2h = gini(length(p)-d + (1:length(q)-d))';
  
  res(d,2) =  norm([p;q] -[conv(hh,g1h)';conv(hh,g2h)'],2);
  [hh, g1h, g2h, resn, iter] = c_f_newton_iter_mod(p',q', ...
                                   hh,g1h,g2h, 1, 100, 1e-15);
  iters(d,4) = iter;
  res(d,7) =  norm([p;q] -[conv(hh,g1h)';conv(hh,g2h)'],2);
  
  d
  [hh, g1h, g2h,  uvres, uvcond] = uvGCDfixedDegree(p,q,d)  
  res(d,4) =  norm([p;q] -[conv(hh,g1h);conv(hh,g2h)],2);
end

save('res_terui5.txt', 'res', '-ascii');
save('iter_terui5.txt', 'iters', '-ascii');