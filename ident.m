% IDENT - LTI system identification by Hankel low-rank approximation
% [sysh, info, wh] = ident(w, m, ell, opt)
%
% W   - T samples, Q variate time series, stored in a TxQ array
%       or N such trajectories, stored in a TxQxN array
%       or N trajectories, stored in cell array with TixQ elements
%       Missing elements are denoted by NaN's.
% M   - number of inputs (P := Q - M outputs)
% ELL - system lag (the order is N := ELL * P)
% OPT - options for the optimization algorithm:
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.WINI = 0 specifies zero initial conditions (default [])
%   OPT.SYS0 - initial approximation: SS object with M inputs, 
%       P outputs, and order N; or Px((ELL + 1) * Q) kernel parameter 
%   OPT.SS = 0 - do not convert the identified model to state space 
%   OPT.DISP - level of display [off|iter|notify|final] (default off)
%   OPT.MAXITER - maximum number of iterations (default 100)
% SYSH - identified system (by default an SS object)
% INFO - information from the optimization solver:
%   INFO.M - misfit ||W - WH||_F^2
%   INFO.ITER - number of iterations
% WH    - optimal approximating time series
function [sysh, info, wh] = ident(w, m, ell, opt)
if ~exist('opt'), opt = []; end
if ~isfield(opt, 'exct'), opt.exct = []; end
if ~iscell(w)
  [T, q, N] = size(w); T = ones(N, 1) * T;
else
  N = length(w); for k = 1:N, [T(k), q] = size(w{k}); end, T = T(:);
end, p = q - m; n = p * ell; 
s.m = (ell + 1) * ones(q, 1); s.n = T - ell; r = (ell + 1) * m + n;
if ~isempty(opt.exct)
  if iscell(opt.exct)
    for i = 1:N, s.w(:, i) = ones(q, 1); s.w(opt.exct{i}) = inf; end
  else 
    s.w = ones(q, 1); s.w(opt.exct) = inf; 
  end
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
  if isfield(s, 'w'), s = rmfield(s, 'w'); end
  if ~iscell(opt.wini) && ~isempty(opt.wini)
    W = ones(T, q, N); W(:, opt.exct, :) = inf;
    s.n = s.n + ell; T = T + ell;
    s.w = [inf * ones(ell, q, N); W]; 
    w   = [opt.wini; w];
  elseif iscell(opt.wini) && ~isempty(cell2mat(opt.wini))
    for k = 1:N
      W = ones(T(k), q); W(:, opt.exct{k}) = inf;
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
if isfield(opt, 'sys0') 
  if isa(opt.sys0, 'lti'),  
    opt.Rini = ss2r(opt.sys0); 
  else, 
    opt.Rini = opt.sys0; 
  end
end
if isfield(s, 'w'),
  if isfield(opt, 'wini'), s.w = w2p(s.w); end
  if all(size(s.w) == [q 1]), s.w = s.w(:, ones(1, N)); end
end
if isfield(opt, 'v'),
  if all(size(opt.v) == [q 1]), opt.v = opt.v(:, ones(1, N)); end
  if ~all(size(opt.v) == [q N]), opt.v = w2p(opt.v); end
end
if isfield(s, 'w') && isfield(opt, 'v')
  try, s.w = s.w .* opt.v; catch
    q = length(s.m); if exist('p'), 
                       np = length(p); 
                     else
                       np = sum(s.m) * length(s.n) + length(s.m) * sum(s.n) ...
                                                   - length(s.m) * length(s.n);
                     end, if ~isfield(s, 'n'), s.n = np - sum(s.m) + 1; end
    N = length(s.n); n  = sum(s.n);
    if length(s.w(:)) == q || all(size(s.w) == [q N])
      % convert q x 1 s.w to q x N
      if isvector(s.w), s.w = s.w(:); s.w = s.w(:, ones(1, N)); end
      % convert q x N s.w to np x 1
      w = [];
      for j = 1:N
        for i = 1:q
          wij = s.w(i, j) * ones(s.m(i) + s.n(j) - 1, 1); w = [w; wij];
        end
      end
      s.w = w;
    end
    w_inf = s.w; s.w = opt.v; 
    q = length(s.m); if exist('p'), 
                       np = length(p); 
                     else
                       np = sum(s.m) * length(s.n) + length(s.m) * sum(s.n) ...
                                                   - length(s.m) * length(s.n);
                     end, if ~isfield(s, 'n'), s.n = np - sum(s.m) + 1; end
    N = length(s.n); n  = sum(s.n);
    if length(s.w(:)) == q || all(size(s.w) == [q N])
      % convert q x 1 s.w to q x N
      if isvector(s.w), s.w = s.w(:); s.w = s.w(:, ones(1, N)); end
      % convert q x N s.w to np x 1
      w = [];
      for j = 1:N
        for i = 1:q
          wij = s.w(i, j) * ones(s.m(i) + s.n(j) - 1, 1); w = [w; wij];
        end
      end
      s.w = w;
    end
    s.w = s.w .* w_inf;
  end
elseif isfield(opt, 'v'), s.w = opt.v; end



[ph, info] = slra(w2p(w), s, r, opt); info.M = info.fmin;
wh = p2w(ph, q, N, T, iscell(w)); 
if isfield(opt, 'ss') && (opt.ss == 0) 
  sysh = info.Rh; 
else 
  sysh = r2ss(info.Rh, m, ell); 
end
if isfield(opt, 'wini')
  if ~iscell(opt.wini) && ~isempty(opt.wini)
    wh = wh(ell + 1:end, :, :);
  elseif iscell(opt.wini) && ~isempty(cell2mat(opt.wini))
    for k = 1:N
      if ~isempty(opt.wini{k})
        wh{k} = wh{k}(ell + 1:end, :);
      end
    end  
  end
end
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
sys = ss(sys); a = sys.a; b = sys.b; c = sys.c; d = sys.d; 
[p, m] = size(d); n = size(a, 1); 
ell1 = n / p + 1; L = ell1; O = c; for t = 2:L, O = [O; O(end - p + 1:end, :) * a]; end, P = null(O')';
if (m > 0)
  F = [d; O(1:(end - p), :) * b]; TT = zeros(ell1 * p, ell1 * m);
  for i = 1:ell1
    TT((i - 1) * p + 1:end, (i - 1) * m + 1: i * m) = F(1:(ell1 + 1 - i) * p, :);
  end
  Q = P * TT; 
else, Q = []; end
R = permute([reshape(Q, p, m, ell1), -reshape(P, p, p, ell1)], [1 3 2]);
R = reshape(R, p, (m + p) * ell1);
function sysh = r2ss(R, m, ell)
[p, tmp] = size(R); ell1 = ell + 1; q = tmp / ell1; n = ell * p;
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
