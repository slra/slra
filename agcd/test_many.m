addpath 'fastgcd';
addpath 'uvgcd';

addpath 'h:/repos/slra';

p1 = [-43.810;-24.300;64.186;23.013;243.81;-195.53;125.80;-266.22;106.68;-185.83;182.73;-16.316];
p2 = [12.517;6.9428;-18.339;-6.5751;-69.659;55.866;-35.944;76.064;-30.481;53.094;-52.209;4.6618];
p3 = [-14.783;14.297;-9.0031;4.3083;15.007;-34.013;83.932;-45.276;2.2157;-59.034;47.507;-4.115];

d = 2;
h130 = lsdivmult({p1,p3}, d, g_ini({p1,p3}, d));
h230 = lsdivmult({p2,p3}, d, g_ini({p2,p3}, d));
h120 = lsdivmult({p1,p2}, d, g_ini({p1,p2}, d));
h1230 = lsdivmult({p1,p2,p3}, d, g_ini({p1,p2,p3}, d));
h1230full = lsdivmult({p1,p2,p3}, d, gini_full_3(p1,p2,p3, d));
[h123mp, sv123] = h_ini_mp({p1,p2,p3}, d,20);
plot(sv123)

hinis = [];
hinis = [h130 h230 h120 h1230 h1230full  h123mp];


figure
[X, Y] = meshgrid(linspace(-3, 3, 100), linspace(-2.5, 2.5, 100));
Z = gcd_nls_evaluate({p1,p2,p3}, [ -ones(1, length(X(:))); X(:)'; Y(:)']);
Z = log(reshape(Z, size(X)));
contour(X,Y,Z);
%surfc(X,Y,Z);
%camorbit(75, 10);
%colorbar;
hold on
plot3([0.96],[-0.08], [0], 'o');

opt.disp = 'iter';
opt.maxiter = 300;
%opt.hini = [-1, 1, -0.4]';

styles = ['x+*sdv'];
names = {'s13','s23','s12','s123', 's123full', 's123mp'}

for i=1:size(hinis,2)
  opt.hini = hinis(:,i);%info23.Rh';
  opt.Psi = [ 0 1 0; 0 0 1; 1 0 0];
  [ph, info] = gcd_nls({p1,p2,p3}, [], d, opt);
  traj = reshape(info.RhK, 3, info.iter+1);
  traj = traj ./ (- ones(3,1) * traj(1,:));
  hold on
  plot(traj(2,:),traj(3,:), ['k' styles(i) '-']);
end  
%plot3(traj(2,:),traj(3,:), log(info.iterinfo(2,:)), 'bo-');




