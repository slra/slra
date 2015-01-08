%% SLRA_DOC - solves the structured low-rank approximation problem
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
%  slra_mex_obj, reg_slra, OptimizationOptions
