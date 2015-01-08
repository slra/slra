%% SLRA - solves the structured low-rank approximation problem
%
%  minimize over ph |p-ph|^2_w subject to rank(S(ph)) <= r
%
%  where  S(ph) is the (m x n)  matrix structure
%         p is a given vector (of length np)
%         r is bound on the rank
%         |p|^2_w = sum(w .* (p.^2)) - the weigthed semi-norm
%             defined by a vector of weights w (of length np)
%             w(i) = Inf  <=>  ph(i) = p(i) (constraint on the approximation)
%
%% Syntax
%   [ph, info] = slra(p, s, r, opt)
%
%% Input
%  p - structure parameter vector
%  s - matrix structure and norm specification (MATLAB structure)
%      The weights are specified by
%         s.w - vector of weights w (default ones(np, 1))
%      The structure can be
%       (a) Mosaic-Hankel-like stucture: S(ph) := Phi * H, where
%             H is a q x N block matrix with Hankel blocks
%             H_ij(ph) = hankel(ph_ij(1:m_i, m_i:(m_i + n_j - 1)))
%             and ph = [ph_11; ... ph_q1; ... ph_1N; ... ph_qN]
%           Defined by:
%             s.m = [m_1 ... m_q], s.n = [n_1 ... n_N] (the sizes of the blocks)
%             s.phi = Phi  (the matrix Phi, by default, identity matrix)
%       (b) A general affine matrix structure S(ph) := S0 + ph(tts)
%           Defined by:
%             s.S0 - the constant m x n matrix (by default, zeros(m,n))
%             s.tts - integer matrix of positions of the elements of ph in S(ph)
%
%  r (optional) - bound on the rank (default is rank reduction by 1)
%  opt (optional) - optimization options (MATLAB structure)
%      opt.Rini - initial approximation (default unstructured LRA)
%                 Rini is a basis for an approximate left kernel of S(p)
%      opt.disp - information about progress of the optimization,
%                 possible values 'off', 'notify', 'iter'
%      opt.solver - solver (default 'c'), can take values:
%          'c' --- efficient C++ solver       (calls SLRA_MEX_OBJ)
%          'm' --- general solver             (calls SLRA_EXT)
%          'r' --- factorization-based solver (calls REG_SLRA)
%
%      ... (additional fields, e.g. opt.psi, depend on the solver being used,
%            see the description of the opt parameter in the solver files help)
%
%% Output
%  ph - approximation structure parameter vector
%  info - information about optimization (MATLAB structure)
%      info.Rh - low-rank certificate Rh (an (m-r) x m matrix)
%                such that Rh * S(ph) = 0
%      info.iter - number of iterations
%      info.time - execution time
%      info.fmin = |p - ph|^2_w - the value of the cost function
%
%
%% See also 
%  slra_mex_obj, reg_slra
function [ph, info] = slra(p, s, r, opt)
if ~exist('opt'), opt = struct; end
if ~isfield(opt, 'solver'), opt.solver = 'c'; end 
if ~isfield(opt, 'disp'), opt.disp = 'off'; end 
if ~isfield(opt, 'tol_m'), opt.tol_m = 1e-6; end 
Im = find(isnan(p));
if ~isempty(Im)
  if ~isfield(s, 'w'), s.w = ones(size(p)); end 
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
  p(Im) = 0; s.w(Im) = 0;
end
if opt.solver == 'c'
  if isfield(s, 'w'), s.w(find(s.w == 0)) = opt.tol_m; end
  opt = rmfield(opt, 'solver');
  obj = slra_mex_obj('new', p, s, r);
  [ph, info] = slra_mex_obj('optimize', obj, opt);
  slra_mex_obj('delete', obj);
elseif opt.solver == 'r'
  if isfield(s, 'w')
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
    opt.w = s.w; 
  end
  if isfield(opt, 'Rini'), opt.P_init = null(opt.Rini); end
  np = length(p); s.tts = s2s(s, np); 
  [ph, info] = reg_slra(p, s, r, opt);
  info.Rh = null(info.P')';
else
  if ~isfield(s, 'w'), s.w = []; end
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
  if ~all(size(s.w) == [length(p(:)) length(p(:))]), s.w = s.w(:); end
  if ~isfield(s,   'phi' ), s.phi = []; end
  if ~isfield(opt, 'Rini'), opt.Rini = []; end
  if ~isfield(opt, 'psi' ), opt.psi = []; end
  np = length(p); warning_state = warning; warning('off');
  [ph, info] = slra_ext(s2s(s, np), p, r, s.w, opt.Rini, s.phi, opt.psi, opt);
  warning('warning_state');
end
function S = s2s(s, np)
q = length(s.m); if exist('p'), 
                   np = length(p); 
                 else
                   np = sum(s.m) * length(s.n) + length(s.m) * sum(s.n) ...
                                               - length(s.m) * length(s.n);
                 end, if ~isfield(s, 'n'), s.n = np - sum(s.m) + 1; end
N = length(s.n); n  = sum(s.n);, p = 1:np;  
if ~isfield(s, 'phi'), s.phi = eye(sum(s.m)); end, [m, mp] = size(s.phi); 
tmp = cumsum([1; s.m(:)]); Imb = tmp(1:end - 1); Ime = tmp(2:end) - 1;
tmp = cumsum([1; s.n(:)]); Inb = tmp(1:end - 1); Ine = tmp(2:end) - 1;
S = zeros(mp, n); ind = 1;
for j = 1:N
  for i = 1:q
    npij = s.m(i) + s.n(j) - 1;
    pij = p(ind:(ind + npij - 1)); ind = ind + npij;
    Hij = hankel(pij(1:s.m(i)), pij(s.m(i):end));
    S(Imb(i):Ime(i), Inb(j):Ine(j)) = Hij;
  end
end
