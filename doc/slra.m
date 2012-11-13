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
% opt.maxiter, opt.tol - stopping criteria
%
% Output arguments:
% ph - approximation structure parameter vector
% info.Rh   - low-rank certificate: Rh * S(ph) = 0
% info.iter - number of iterations
% info.time - execution time
% info.fmin = norm(w .* (p - ph)) ^ 2
%
% Note: it is required that length(p) > n * (m - r).

function [ph, info] = slra(p, s, r, varargin)
non_negative = @(x) isnumeric(x) && all(all(x >= 0));
integer = @(x) isnumeric(x) && ~isempty(x) && all(all(~mod(x, 1)));
ipr = inputParser;
ipr.addRequired('p', @(p) isnumeric(p) && ~isempty(p)); 
ipr.addRequired('s', @(s) isfield(s, 'm'));
ipr.addRequired('r', @(r) integer(r) && r >= 0); 
ipr.parse(p, s, r); 
if ~isfield(s, 'phi'), s.phi = eye(sum(s.m)); end
if ~isfield(s, 'n'), s.n = (length(p) - sum(s.m)) + 1; end
if ~isfield(s, 'w'), s.w = ones(size(p)); end
ip = inputParser; ip.KeepUnmatched = true;
m = size(s.phi, 1); ip.addParamValue('psi', [], @(psi) isnumeric(psi)); 
ip.addParamValue('Rini', [], @(Rini) isnumeric(Rini) && ...
                                     all(size(Rini) == [m - r m]));
ip.addParamValue('solver', 'c', @(solver) solver == 'c' || solver == 'm'); 
ip.addParamValue('method', 'll', ...
       @(method) ischar(method) && length(method) <= 3) % && ...
       %any(strcmpi(method(1), ['reg' combine('lqn', {'ls', 'b2pf', 'n2r'})])));
ip.addParamValue('disp', 'off', ...
       @(disp) any(strcmpi(disp, {'iter', 'notify', 'off'})));
ip.addParamValue('maxiter', 100, @(maxiter) integer(maxiter) && maxiter >= 0);
ip.addParamValue('tol', 1e-5, non_negative);
ip.addParamValue('epsrel', 1e-5, non_negative);
ip.addParamValue('epsabs', 1e-5, non_negative);
ip.addParamValue('epsgrad', 1e-5, non_negative);
ip.addParamValue('step', 0.001, non_negative);
ip.addParamValue('reggamma', 0.001, non_negative); 
ip.parse(varargin{:}); opt = ip.Results;
[m, mp] = size(s.phi); q = length(s.m); N = length(s.n); n = sum(s.n); 
s2np = @(s) sum(s.m) * length(s.n) + length(s.m) * sum(s.n) ...
                                   - length(s.m) * length(s.n);, np = s2np(s); % = N * mp + q * n - q * N;
if opt.solver == 'c', opt.psi = eye(m); else, eye((m - r) * m); end
if sum(size(s.w)) == np + 1, Im = unique([find(isnan(p)) find(s.w == 0)]); w(Im) = 0; p(Im) = NaN; else, Im = []; end
if opt.solver == 'c' && ~isempty(Im), s.w(Im) = 1e-6; p(Im) = 0; end  
if opt.solver == 'c'
  [ph, info] = slra_mex(p, s, r, opt); 
else
  [ph, info] = slra_ext(s2s(s), p, r, diag(s.w), opt.Rini, s.phi, opt.psi, opt); 
end
function c = combine(f, s)
c = {};
for i = 1:length(f)
    c = [c f(i)];
    for j = 1:length(s{i})
        c = [c [f(i) s{i}(j)]];
    end
end
