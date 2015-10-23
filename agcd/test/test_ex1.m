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

res = zeros(n,5);
iters = zeros(n,5);
res(:,1) = (1:n)';
iters(:,1) = (1:n)';

opt.disp = 'iter';
opt.epsgrad = 1e-20;

methods = {'gcd_nothing', 'gcd_nls', 'gcd_cofe', 'gcd_syl', 'gcd_uvgcd', 'gcd_fastgcd'};
names = { '', 'd', 'LRA', 'VP$_{h}$', 'VP$_{g}$', 'VP$_{S}$', 'UVGCD', 'FASTGCD'};

%% Run methods
for d=1:n
  d
  opt.gini = g_ini({p,q}, d);
  opt.hini = lsdivmult({p,q}, d, opt.gini);

  for j=1:length(methods)
    eval(['[ph, info] = ' methods{j} '({p,q}, [], d, opt);']) 
    res(d, j+1) = norm([p;q] -cell2mat(ph),2);
    iters(d, j+1) = info.iter;  
  end    
end  
%res(:,3) = [5.17e-1;6.95e-4;1.97e-5;2.89e-6;5.28e-5;...
%            2.15e-3;8.34e-2;2.04e0;4.70e1;7.73e2] / sccoef;

fid = fopen('res_ex1_matlab.txt', 'w');
fprintf(fid, '\t%s', names{:});
fprintf(fid, '\n');
fclose(fid);
save('res_ex1_matlab.txt', 'res', '-ascii','-append');
fid = fopen('iters_ex1_matlab.txt', 'w');
fprintf(fid, '\t%s', names{:});
fprintf(fid, '\n');
fclose(fid);
save('iters_ex1_matlab.txt', 'iters', '-ascii','-append');
