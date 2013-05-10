import manopt.solvers.trustregions.*;
import manopt.manifolds.sphere.*;
import manopt.tools.*;

T = 400;
t = 1:T;
p0 = 2*cos((2*pi*t)/3) + cos((2*pi*t)/7) + cos((2*pi*t)/10);
p =  p0 + randn(1, T) * 0.02 * i;

r = 6;
s.m = r+1;

obj = slra_mex_obj('new', p,s,r);

% Create the problem structure.
manifold = grassmannfactory(s.m) s.m-r);
problem.M = manifold;
    
% Define the problem cost function and its gradient.
problem.cost = @(x) slra_mex_obj('func', obj, x');
problem.grad = @(x) manifold.proj(x, slra_mex_obj('grad', obj, x')');
    
% Numerically check gradient consistency.
checkgradient(problem);
    
slra_mex_obj('delete', obj);    
    

