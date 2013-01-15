% IDENT - Identification by structured low-rank approximation
%
% [sysh, info, wh, xini] = ident(w, m, ell, opt)
%
% Inputs:
% W   - set of trajectories, stored in an arry with dimensions
%                 #samples x #variables x #time series        
%       in case of equal number of samples or a cell array with  
%       #time series entries of dimension #samples x #variables
% M   - input dimension
% ELL - system lag
% OPT - options for the optimization algorithm:
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.WINI = 0 specifies zero initial conditions (default [])
%   OPT.SYS0 - initial approximation: an SS system with M inputs, 
%              P := size(W, 2) - M outputs, and order N : = ELL * P
%   OPT.DISP - level of display [off | iter | notify | final] (default notify)
%   OPT.MAXITER - maximum number of iterations (default 100)
%   OPT.TOL  - convergence tolerance (default 1e-4)
%
% Outputs:
% SYSH - input/state/output representation of the identified system
% INFO - information from the optimization solver:
%   INFO.M - misfit ||W - WH||_F,
%   INFO.TIME - execution time for the SLRA solver,
%   INFO.ITER - number of iterations
% W    - optimal approximating time series
% XINI - initial conditions under which WH are obtained
function [sysh, info, wh, xini] = ident(w, m, ell, varargin)
ipr = inputParser;
ipr.addRequired('w', @(w) (isnumeric(w) && ~isempty(w)) || iscell(w)); 
ipr.addRequired('m', @(m) m >= 0);
ipr.addRequired('ell', @(ell) ell >= 0); 
ipr.parse(w, m, ell); 
if ~iscell(w)
  [T, q, N] = size(w); T = ones(N, 1) * T;
else
  N = length(w); for k = 1:N, [T(k), q] = size(w{k}); end, T = T(:);
end
p = q - m; n = p * ell; ell1 = ell + 1;
ip = inputParser; ip.KeepUnmatched = true;
ip.addParamValue('exct', [], @(exct) isempty(exct) || ...
                                     (all(exct >= 1 & exct <= q)))
ip.addParamValue('wini', [])
ip.addParamValue('sys0', [], @(sys0) isempty(sys0) || isa(sys0, 'lti'));
ip.parse(varargin{:}); opt = ip.Results;
ip.addParamValue('disp', 'off');
ip.addParamValue('method', 'll');
ip.addParamValue('maxiter', 100);
ip.addParamValue('tol', 1e-5);
ip.parse(varargin{:}); opt = ip.Results;
if isa(opt.sys0, 'lti'), [a, b, c, d] = ssdata(ss(opt.sys0)); L = ell1; O = c; for t = 2:L, O = [O; O(end - size(c, 1) + 1:end, :) * a]; end
                         R = null(O')';
                         if (m > 0)
                           F = [d; O(1:(end - p), :) * b];
                           TT = zeros(ell1 * p, ell1 * m);
                           for i = 1:ell1
                             TT((i - 1) * p + 1:end, (i - 1) * m + 1: i * m) = F(1:(ell1 + 1 - i) * p, :);
                           end
                           Q = R * TT;
                           R3 = [reshape(Q, p, m, ell1), -reshape(R, p, p, ell1)];
                           R  = reshape(R3, p, ell1 * q);
                         end , opt.Rini = R; end
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
if ~iscell(w)
  par = vec([opt.wini; w]);
else
  par = []; for k = 1:N, par = [par; vec([opt.wini{k}; w{k}])]; end
end
s.m = ell1 * ones(q, 1); s.n = T - ell; r = ell1 * m + n;
L = q * ell1; phi = []; 
for i = 1:ell1, phi = [phi, i:ell1:L]; end
s.phi = eye(L); s.phi = s.phi(phi, :);
if ~isempty(opt.exct), s.w = ones(q, 1); s.w(opt.exct) = inf; end
if ~iscell(opt.wini) && ~isempty(opt.wini)
  W = ones(T, q, N); W(:, opt.exct, :) = inf;
  s.n = s.n + ell; T = T + ell;
  s.w = vec([inf * ones(ell, q, N); W]);
elseif iscell(opt.wini) && ~isempty(cell2mat(opt.wini))
  s.w = []; 
  for k = 1:N
    W = ones(T(k), q); W(:, opt.exct) = inf;
    if ~isempty(opt.wini{k})
      s.n(k) = s.n(k) + ell; T(k) = T(k) + ell;
      s.w = [s.w; vec([inf * ones(ell, q); W])];
    else
      s.w = [s.w; W(:)];
    end
  end
end
Im = find(isnan(par));
if ~isempty(Im)
  if ~isfield(s, 'w'), s.w = ones(size(par));
  elseif all(size(s.w) == [q 1]) 
    if iscell(w)
      s.w = []; 
      for k = 1:N
        W = ones(T(k), q); W(:, opt.exct) = inf; s.w = [s.w; W(:)];
      end
    else
      W = ones(T, q, N); W(:, opt.exct, :) = inf; s.w = vec(W);
    end
  end
  s.w(Im) = 1e-5; par(Im) = 0;
end
[ph, info] = slra(par, s, r, opt); info.M = info.fmin;
if ~iscell(w)
  wh = reshape(ph, T(1), q, N); 
  if ~isempty(opt.wini), wh = wh(ell1:end, :, :); end
else
  for k = 1:N, 
    wh{k} = reshape(ph(1:(q * T(k))), T(k), q); ph(1:(q * T(k))) = [];
    if ~isempty(opt.wini{k}), wh{k} = wh{k}(ell1:end, :); end
  end
end
info.Rh = - info.Rh / info.Rh(:, end - p + 1:end);
if m > 0 
  R3 = reshape(info.Rh, p, q, ell1);
  Q3 = R3(:, 1:m, :); P3 = - R3(:, m + 1:q, :);
  a = zeros(n); b = zeros(n, m); c = [];
  if n > 0
    a(p + 1:end, 1:n - p) = eye(n - p); 
    c = [zeros(p, n - p) eye(p)];
  end
  d = Q3(:, :, ell1); ind_j = (n - p + 1):n;
  for i = 1:ell
    ind_i = ((i - 1) * p + 1):(i * p); P3i = P3(:, :, i);
    a(ind_i, ind_j) = - P3i; 
    b(ind_i, :) = Q3(:, :, i) - P3i * d;
  end
  sysh = ss(a, b, c, d, 1);
else 
  O = null(info.Rh); 
  a = O(1:end - p, :) \ O(p + 1:end, :); c = O(1:p, :);
  sysh = ss(a, [], c, [], 1); 
end
if nargout > 3, xini = inistate(wh, sysh); end
function xini = inistate(w, sys)
[a, b, c, d] = ssdata(sys); [p, m] = size(d); n = size(a, 1); 
if ~iscell(w)
  [T, q, N] = size(w); T = ones(N, 1) * T;
else
  N = length(w); for k = 1:N, [T(k), q] = size(w{k}); end, T = T(:);
end, L = max(T); 
O = c; for t = 2:L, O = [O; O(end - size(c, 1) + 1:end, :) * a]; end, sys.Ts = -1; xini = zeros(n, N);  
for k = 1:N
  if ~iscell(w)
    uk = w(:, 1:m, k); yk = w(:, (m + 1):end, k); Tk = T(1);
  else  
    uk = w{k}(:, 1:m); yk = w{k}(:, (m + 1):end); Tk = T(k);
  end
  if m > 0, y0k = (yk - lsim(sys, uk, 0:(Tk - 1)))'; else, y0k = yk'; end
  xini(:, k) = O(1:(Tk * p), :) \ y0k(:);
end
function a = vec(A), a = A(:);
