% MISFIT - orthogonal distance from W to an LTI system SYS
% M = minimum over Wh  ||W - Wh||_F^2 s.t. Wh is a trajectory of SYS
% [M, wh] = misfit(w, sys, opt)
% 
% W   - time series or a set of time series (see "help ident")
% SYS - LTI system: SS object with M inputs, P := Q - M outputs and 
%       order N := ELL * P, or Px((ELL + 1) * Q) kernel parameter 
% OPT - options 
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.WINI = 0 specifies zero initial conditions (default [])
% M    - misfit ||W - WH||_F^2
% WH   - optimal approximating time series
function [M, wh] = misfit(w, sysh, opt)
if isa(sysh, 'ss')
  [p, m] = size(sysh); n = size(sysh, 'order'); ell = n / p;
else 
  [p, ell1q] = size(sysh); 
  if ~iscell(w)
    [T, q, N] = size(w); T = ones(N, 1) * T;
  else
    N = length(w); for k = 1:N, [T(k), q] = size(w{k}); end, T = T(:);
  end
  ell = ell1q / q - 1; m = q - p;
end  
opt.sys0 = sysh; opt.maxiter = 0; opt.disp = 'off';
[~, info, wh] = ident(w, m, ell, opt); M = info.M;
