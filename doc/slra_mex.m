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
function varargout = slra_mex(p, s, r, opt)
  obj = slra_mex_obj('new', p, s, r);
  [varargout{1:nargout}] = slra_mex_obj('optimize', obj, opt);
  slra_mex_obj('delete', obj);
end
