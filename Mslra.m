function [ph, info] = slra_misfit(p, s, R, psi)
opt.maxiter = 0; opt.Rini = R; 
if nargin == 4, opt.psi = psi; end
r = size(R, 2) - size(R, 1); 
[ph, info] = slra(p, s, r, opt);
