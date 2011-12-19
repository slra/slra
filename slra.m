% SLRA - solves the structured low rank approximation problem 
% 
% minimize over Ph (P - Ph)' * (P - Ph) subject to rank(S(Ph)) <= r,
%
% where the matrix S(Ph) is of the form [S1 ... Sq], with Si being 
% (block) Toeplitz, Hankel structured, unstructured, or constant
%
% [xh, info, vxh, ph] = slra(p, s, r, x0, opt, phi)
%
% Input arguments:
%
% P - structure parameter vector
% S - structure specification: a matrix with q rows and 1-4 columns:
%   S(i, 1) - number of subblocks in a block row of Si 
%             (use S(i, 1) = 1 for unstructured block Si)
%   S(i, 2) - number of columns in a subblock (default 1)
%   S(i, 3) = 1 if the block is exact, or 0 otherwise (default 0)
%   S(i, 4) = 1 if the block is Toeplitz, 0 if it is Hankel (default 0)
% For block-Toeplitz/Hankel structured block Si with K x nu(i) size blocks
% S is a structure with fields K and A, where S.A is the matrix described 
% above and S.K = K (the same for all block-Toeplitz/Hankel blocks S1)
% R - rank (default is rank reduction by 1, i.e., R = sum(S(:,2)) - 1) 
% X0  - initial approximation (default unstructured low rank approximation)
%       PHI * [X0; -I] is a basis for an approximate kernel of S(P) 
% OPT - optimization parameters, structure with fields: 
%   OPT.METHOD - two letter string, specifying the optimization method
%       first letter options: 'l' - Levenberg-Marquardt  (default)
%                             'q' - Quasi-Newton, 'n' - Nelder-Mead
%       second letter options for Levenberg-Marquardt:
%                             'l' - lmder (default), 's' - lmsder
%       second letter options for Quasi-Newton:
%       'b' - bfgs, '2' - bfgs2, 'f' - conjugate_fr, 'p' - conjugate_pr
%       second letter options for Nelder-Mead:
%       '?' - nmsimplex, '?' - nmsimplex2, '?' - nmsimplex2rand
%   OPT.MAXITER - maximum number of iterations 
%   OPT.EPSREL, OPT.EPSABS, OPT.EPSGRAD - convergence tolerances
%   OPT.REGGAMMA - regularization parameter for gamma, absolute
%   OPT.DISP - information about progress of the optimization 
%              options: 'iter', 'notify', 'final', 'off' (default 'notify')
% PHI - transformation matrix (default I)
%
% Output arguments:
%
% XH   - low rank certificate: S(Ph) * PHI * [XH; -I] = 0
% INFO - convergence information: 
%   INFO.ITER - number of iterations
%   INFO.TIME - execution time
%   INFO.FMIN = (P - Ph)' * (P - Ph)
% VXH  - asymptotic covariance matrix of vec(XH)
%
% Note: it is required that length(P) > size(S(P), 1) * (size(S(P), 2) - R)
