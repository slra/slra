% IDENT - LTI system identification by structured low-rank approximation
% [sysh, info, wh, xini] = ident(w, m, ell, opt)
%
% W   - T samples, Q variate time series, stored in a TxQ array
%       or N such trajectories, stored in a TxQxN array
%       or N trajectories, stored in cell array with TixQ elements
%       Missing elements are denoted by NaN's.
% M   - input dimension
% ELL - system lag (order = ELL * (Q - M))
% OPT - options for the optimization algorithm:
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.WINI = 0 specifies zero initial conditions (default [])
%   OPT.SYS0 - initial approximation: an SS system with M inputs, 
%              P := Q - M outputs, and order N : = ELL * P
%   OPT.DISP - level of display [off|iter|notify|final] (default off)
%   OPT.MAXITER - maximum number of iterations (default 100)
% SYSH - input/state/output representation of the identified system
% INFO - information from the optimization solver:
%   INFO.M - misfit ||W - WH||_F^2
%   INFO.ITER - number of iterations
% WH    - optimal approximating time series
% XINI - initial condition under which WH is obtained by SYSH
function [sysh, info, wh, xini] = ident(w, m, ell, opt)
if ~exist('opt'), opt = []; end
if ~iscell(w)
  [T, q, N] = size(w); T = ones(N, 1) * T;
else
  N = length(w); for k = 1:N, [T(k), q] = size(w{k}); end, T = T(:);
end, p = q - m; n = p * ell; 
s.m = (ell + 1) * ones(q, 1); s.n = T - ell; r = (ell + 1) * m + n;
if isfield(opt, 'exct') && ~isempty(opt.exct) 
  s.w = ones(q, 1); s.w(opt.exct) = inf; 
end
if isfield(opt, 'wini')
  if iscell(w)
    if isempty(opt.wini) 
      opt.wini = cell(1, N); 
    elseif ~iscell(opt.wini) & (opt.wini == 0)
      for k = 1:N, wini{k} = zeros(ell, q); end, opt.wini = wini;
    else  
      for k = 1:N 
        if opt.wini{k} == 0, opt.wini{k} = zeros(ell, q); end
      end 
    end
  elseif opt.wini == 0
    opt.wini = zeros(ell, q, N);
  end
  if ~iscell(opt.wini) && ~isempty(opt.wini)
    W = ones(T, q, N); W(:, opt.exct, :) = inf;
    s.n = s.n + ell; T = T + ell;
    s.w = [inf * ones(ell, q, N); W]; 
    w   = [opt.wini; w];
  elseif iscell(opt.wini) && ~isempty(cell2mat(opt.wini))
    for k = 1:N
      W = ones(T(k), q); W(:, opt.exct) = inf;
      if ~isempty(opt.wini{k})
        s.n(k) = s.n(k) + ell; T(k) = T(k) + ell;
        s.w{k} = [inf * ones(ell, q); W]; 
        w{k}   = [opt.wini{k}; w{k}];
      else
        s.w{k} = W;
      end
    end  
  end
end
if isfield(opt, 'sys0') && isa(opt.sys0, 'lti'), opt.Rini = ss2r(opt.sys0); end
par = w2p(w); if isfield(s, 'w') && length(s.w) == length(par), s.w = w2p(s.w); end

save example

[ph, info] = slra(par, s, r, opt); info.M = info.fmin;
wh = p2w(ph, q, N, T, iscell(w)); sysh = r2ss(info.Rh, m, ell); 
if nargout > 3, xini = inistate(wh, sysh); end
function p = w2p(w)
if ~iscell(w), p = w(:); else
  p = []; for k = 1:length(w), p = [p; w2p(w{k})]; end
end
function w = p2w(p, q, N, T, c)
if ~c 
  if N == 1, w = reshape(p, T, q, N); else
    for k = 1:N, w(:, :, k) = p2w(p(1:q * T(k)), q, 1, T(k), 0); p(1:q * T(k)) = []; end
  end
else
  for k = 1:N, w{k} = p2w(p(1:q * T(k)), q, 1, T(k), 0); p(1:q * T(k)) = []; end
end
function R = ss2r(sys)
[a, b, c, d] = ssdata(ss(sys)); [p, m] = size(d); n = size(a, 1); 
ell1 = n / p + 1; L = ell1; O = c; for t = 2:L, O = [O; O(end - p + 1:end, :) * a]; end, P = null(O')';
if (m > 0)
  F = [d; O(1:(end - p), :) * b];
  TT = zeros(ell1 * p, ell1 * m);
  for i = 1:ell1
    TT((i - 1) * p + 1:end, (i - 1) * m + 1: i * m) = F(1:(ell1 + 1 - i) * p, :);
  end
  Q = P * TT; 
else, Q = []; end
R = permute([reshape(Q, p, m, ell1), -reshape(P, p, p, ell1)], [1 3 2]);
R = reshape(R, p, (m + p) * ell1);
function sysh = r2ss(R, m, ell)
[p, tmp] = size(R); ell1 = ell + 1; q = tmp / ell1; n = ell * p;
if m > 0 
  R = permute(reshape(R, p, ell1, q), [1 3 2]);
  Q = R(:, 1:m, :); P = - R(:, m + 1:q, :); inv_Pl = pinv(P(:, :, ell + 1));
  a = zeros(n); b = zeros(n, m); c = [];
  if n > 0
    a(p + 1:end, 1:n - p) = eye(n - p); 
    c = [zeros(p, n - p) eye(p)];
  end
  d = inv_Pl * Q(:, :, ell1); ind_j = (n - p + 1):n;
  for i = 1:ell
    ind_i = ((i - 1) * p + 1):(i * p); Pi = P(:, :, i);
    a(ind_i, ind_j) = - inv_Pl * Pi; 
    b(ind_i, :) = inv_Pl * (Q(:, :, i) - Pi * d);
  end
  sysh = ss(a, b, c, d, 1);
else 
  O = null(R); a = O(1:end - p, :) \ O(p + 1:end, :); c = O(1:p, :);
  sysh = ss(a, [], c, [], 1); 
end
function xini = inistate(w, sys, use_all_data)
[a, b, c, d] = ssdata(sys); [p, m] = size(d); n = size(a, 1); 
if ~iscell(w)
  [T, q, N] = size(w); T = ones(N, 1) * T;
else
  N = length(w); for k = 1:N, [T(k), q] = size(w{k}); end, T = T(:);
end 
if ~exist('use_all_data') || use_all_data ~= 1, T = n * ones(1, N); end
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
