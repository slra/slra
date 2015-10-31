clear all;
addpath '..';
addpath 'uvgcd';
addpath 'fastgcd';

n = 10;

%% Generate polynomials
p = 1; q = 1;
for j=1:n
  x_j = (-1)^j * j/2;
  p = conv(p,[-x_j; 1]);
  q = conv(q,[-x_j+ 10^(-j); 1]);
end

sccoef = norm(p,2);
p = p./ norm(p,2); %sccoef;
q = q./ norm(q,2); %sccoef;


%opt.disp = 'iter';
opt.epsgrad = 1e-20;

methods = {'gcd_nothing', 'gcd_nls', 'gcd_cofe', 'gcd_syl', 'gcd_uvgcd', 'gcd_fastgcd'};
names = {'LRA', 'VP$_{h}$', 'VP$_{g}$', 'VP$_{S}$', 'UVGCD', 'FASTGCD'};

res = zeros(n,length(methods));
iters = zeros(n,length(methods));


%% Run methods
for d=1:n
  d
  opt.gini = g_ini({p,q}, d);
  opt.hini = lsdivmult({p,q}, d, opt.gini);

  for j=1:length(methods)
    eval(['[ph, info] = ' methods{j} '({p,q}, [], d, opt);']) 
    res(d, j) = norm([p;q] -cell2mat(ph),2);
    iters(d, j) = info.iter;  
  end    
end  

conds = zeros(n,3);

for d=1:n
  ell = n-d;  
  opt.gini = g_ini({p,q}, d);
  opt.hini = lsdivmult({p,q}, d, opt.gini);
  [ph, info] = gcd_uvgcd({p,q}, [], d, opt);
  
  gh1 = info.bfgh{1}; gh2 = info.bfgh{2};
  conds(d,1) = cond(multmat(info.hh, ell))^2;
  conds(d,2) = cond([multmat(gh1, d); multmat(gh2, d)])^2;
  
  A = multmat(gh1,n+ell+1).'; B=  multmat(gh2,n+ell+1).';
  A
  A(:,ell+1:ell+2+n)
  conds(d,3) = cond([A(:,ell+1:ell+2+n) B(:,ell+1:ell+2+n)])^2;
end

%res(:,3) = [5.17e-1;6.95e-4;1.97e-5;2.89e-6;5.28e-5;...
%            2.15e-3;8.34e-2;2.04e0;4.70e1;7.73e2] / sccoef;

save_text_table('res_ex1_matlab.txt', [{'d'} names], [(1:n)' res]);
save_text_table('iters_ex1_matlab.txt', [{'d'} names], [(1:n)' iters]);
save_text_table('conds_ex1_matlab.txt', {'d','VP$_{h}$', 'VP$_{g}$', 'VP$_{S}$'}, [(1:n)' conds]);

