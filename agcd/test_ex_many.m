p1 = [-43.810;-24.300;64.186;23.013;243.81;-195.53;125.80;-266.22;106.68;-185.83;182.73;-16.316];
p2 = [12.517;6.9428;-18.339;-6.5751;-69.659;55.866;-35.944;76.064;-30.481;53.094;-52.209;4.6618];
p3 = [-14.783;14.297;-9.0031;4.3083;15.007;-34.013;83.932;-45.276;2.2157;-59.034;47.507;-4.115];

d = 2;
h130 = lsdivmult({p1,p3}, d, g_ini({p1,p3}, d));
h230 = lsdivmult({p2,p3}, d, g_ini({p2,p3}, d));
h120 = lsdivmult({p1,p2}, d, g_ini({p1,p2}, d));
h1230 = lsdivmult({p1,p2,p3}, d, g_ini({p1,p2,p3}, d));
h1230full = lsdivmult({p1,p2,p3}, d, gini_full_3(p1,p2,p3, d));
[h123mp, sv123] = h_ini_mp({p1,p2,p3}, d,10);
href = [11.64469379842480; -11.28371806974011; 1];
hinis = [h130 h230 h120 h1230 h1230full  h123mp href];
hinis = hinis ./ -(ones(d+1,1) * hinis(1,:));

opt.disp = 'off';
opt.maxiter = 300;
opt.psi = [ 0 1 0; 0 0 1; 1 0 0];

res_table = zeros(2 * d + 2, size(hinis, 2));
res_table(1:d, :) = hinis(2:d+1, :);
infos = cell(1, size(hinis, 2));
for i=1:length(infos)
  opt.hini = hinis(:,i);%info23.Rh';
  opt.maxiter = 0;
  
  [ph, info] = gcd_nls({p1,p2,p3}, [], d, opt);
  opt.maxiter = 300;
  [ph, infos{i}] = gcd_nls({p1,p2,p3}, [], d, opt);
  res_table(d+1, i) = info.fmin;
  res_table(d+2, i) = norm(hinis(:,i) - hinis(:,end), 2);
  res_table(d+2+(1:d), i) = infos{i}.Rh(2:d+1);
  res_table(2*d+3, i) = infos{i}.fmin;
  res_table(2*d+4, i) = infos{i}.iter;
end
save('res_ex3.txt', 'res_table', '-ascii');


figure
[X, Y] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));
Z = gcd_nls_evaluate({p1,p2,p3}, [ -ones(1, length(X(:))); X(:)'; Y(:)']);
Z = log(reshape(Z, size(X)));
surfc(X,Y,Z);
camorbit(20, -10);
