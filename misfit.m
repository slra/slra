% MISFIT - orthogonal distance from W to an LTI system SYS
% M = minimum over Wh  ||W - Wh||_F^2 s.t. Wh is a trajectory of SYS
% [M, wh, xini] = misfit(w, sys, opt)
% 
% W   - time series or a set of time series (see "help ident")
% SYS - state-space system with M inputs, P := size(W, 2) - M outputs
%        and order N := L * P, where L is an integer
% OPT - options 
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.WINI = 0 specifies zero initial conditions (default [])
% M    - misfit ||W - WH||_F^2
% WH   - optimal approximating time series
% XINI - initial condition under which WH is obtained by SYSH
function [M, wh, xini] = misfit(w, sysh, opt)
[p, m] = size(sysh); n = size(sysh, 'order'); ell = n / p;
opt.sys0 = sysh; opt.maxiter = 0; opt.disp = 'off';
[~, info, wh] = ident(w, m, ell, opt); M = info.M;
if nargout > 2, xini = inistate(wh, sysh); end
function xini = inistate(w, sys, use_all_data)
a = sys.a; c = sys.c; [p, m] = size(sys.d); n = size(a, 1); 
if ~iscell(w)
  [T, q, N] = size(w); T = ones(N, 1) * T;
else
  N = length(w); for k = 1:N, [T(k), q] = size(w{k}); end, T = T(:);
end 
if ~exist('use_all_data') || use_all_data ~= 1, T = max(n, 2) * ones(1, N); end
L = max(T); O = c; for t = 2:L, O = [O; O(end - p + 1:end, :) * a]; end
sys.Ts = -1; xini = zeros(n, N);  
for k = 1:N
  if ~iscell(w)
    uk = w(1:T(k), 1:m, k); yk = w(1:L, (m + 1):end, k);
  else  
    uk = w{k}(1:T(k), 1:m); yk = w{k}(1:L, (m + 1):end);
  end
  if m > 0, y0k = (yk - lsim(sys, uk, 0:(T(k) - 1)))'; else, y0k = yk'; end
  xini(:, k) = O(1:(T(k) * p), :) \ y0k(:);
end
