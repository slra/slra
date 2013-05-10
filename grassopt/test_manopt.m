     
% Generate the problem data.
T = 100;
p = 1:T;
p = p + 0.05 * randn(1, T)
r = 2;
s = struct('m', r+1, 'n', T-r);

[ph, info] = slra_grass(p, s, r, struct('maxiter', 100, 'disp', 'iter'));
     
% Create the problem structure.
%manifold = grassmannfactory(s.m, s.m-r,1);
%problem.M = manifold;
     
% Define the problem cost function and its gradient.
%problem.cost = @(x) slra_mex_obj('func', obj, x');
%problem.grad = @(x) manifold.proj(x, slra_mex_obj('grad', obj, x')');
     
% Numerically check gradient consistency.
%checkgradient(problem);
     
% Solve.
%[x xcost info] = trustregions(problem);
     
% Display some statistics.
%figure;
%loglog([info.time], [info.cost], '.-');

%slra_mex_obj('delete', obj);
