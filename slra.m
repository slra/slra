% SLRA - solves the structured low-rank approximation problem 
% minimize over ph norm(w .* (p - ph)) ^ 2 subject to rank(S(ph)) <= r
% where S(ph) = Phi * H, with H a q x N block matrix with Hankel blocks 
% H_ij(ph) = hankel(ph_ij(1:m_i, m_i:(m_i + n_j - 1))) 
% and ph = [ph_11; ... ph_q1; ... ph_1N; ... ph_qN].
%
% [ph, info] = slra(p, s, r, opt)
%
% Input arguments:
% p - structure parameter vector
% s - problem structure specification: 
%     s.m = [m_1 ... m_q], s.n = [n_1 ... n_N], 
%     s.phi = Phi (default identity matrix)
%     s.w = w (default ones(np, 1))
%     (w(i) == inf <=> ph(i) == p(i), w(i) == 0 <=> p(i) not used)
% r - rank (default is rank reduction by 1)
% opt.psi - kernel basis constraint (default identity matrix)
% opt.Rini - initial approximation (default unstructured LRA)
%    Rini is a basis for an approximate left kernel of S(p) 
% opt.disp - information about progress of the optimization 
% opt.solver - solver: 'c' --- efficient, 'm' --- general (default 'c')
% opt.method - optimization method
%
% Output arguments:
% ph - approximation structure parameter vector
% info.Rh   - low-rank certificate: Rh * S(ph) = 0
% info.iter - number of iterations
% info.time - execution time
% info.fmin = norm(w .* (p - ph)) ^ 2
%
% Note: it is required that length(p) > n * (m - r).
function varargout = slra(p, s, r, opt)    

tol_missing = 1e-6; % tolerance for missing value weights in the C solver

% default solver
if ~exist('opt'), opt = []; end
if isfield(opt, 'solver') solver = opt.solver; else solver = 'c'; end

% convert NaNs in p to 0s and create or modify s.w accordingly
Im = find(isnan(p));
if ~isempty(Im)
  p(Im) = 0;
  if ~isfield(s, 'w'), 
    s.w = ones(size(p)); 
  elseif length(s.w) ~= length(p) % convert s.w to np x 1 vector
    % define constants
    q = length(s.m); np = length(p);
    if ~isfield(s, 'n'), s.n = (np - sum(s.m)) + 1; end % default s.n
    N = length(s.n); n = sum(s.n); 
    % convert q x 1 s.w to q x N
    if isvector(s.w), s.w = s.w(:); s.w = s.w(:, ones(1, N)); end
    % convert q x N s.w to np x 1
    w = [];
    for j = 1:N
      for i = 1:q 
        wij = s.w(i, j) * ones(s.m(i), s.n(j)); w = [w; wij(:)]; 
      end 
    end
    s.w = w;
  end
  s.w(Im) = 0;
end
if solver == 'c' && isfield(s, 'w'), s.w(find(s.w == 0)) = tol_missing; end

% call the solver
if solver == 'c'
  obj = slra_mex_obj('new', p, s, r);
  [varargout{1:nargout}] = slra_mex_obj('optimize', obj, opt);
  slra_mex_obj('delete', obj);
else
  if ~isfield(s, 'w'),      s.w = [];      end
  if ~isfield(s, 'phi'),    s.phi = [];    end
  if ~isfield(s, 'w'),      s.w = [];      end
  if ~isfield(opt, 'Rini'), opt.Rini = []; end
  if ~isfield(opt, 'psi'),  opt.psi = [];  end
  dir_str = fileparts(which('slra')); addpath([dir_str '/doc']);
  [varargout{1:nargout}] = slra_ext(s2s(s), p, r, s.w, opt.Rini, s.phi, opt.psi, opt); 
end
