% SLRA - solves the structured low rank approximation problem 
% 
% minimize over Ph ||P - Ph||_2 subject to rank(S(Ph)) <= r,
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
%
% R - rank (default is rank reduction by 1, i.e., R = sum(S(:,1).*S(:,2))-1) 
%
% X0  - initial approximation (default unstructured low rank approximation)
%       PHI * [X0; -I] is a basis for an approximate kernel of S(P) 
%
% OPT - optimization parameters, structure with fields: 
%  Most widely used: 
%    OPT.DISP - information about progress of the optimization 
%               options: 'iter', 'notify', 'final', 'off' (default 'notify')
%    OPT.MAXITER  - maximum number of iterations, 
%    OPT.METHOD   - optimization method (consult GSL library manual): 
%        'l' - Levenberg-Marquardt, gsl_multifit_fdf_solver_...
%        'll'  - ...lmder (default)
%        'ls'  - ...lmsder
%        'q' - minimization with derivatives, gsl_multimin_fdfminimizer_...
%        'qb'  - ...bfgs (default)
%        'q2'  - ...bfgs2
%        'qp'  - ...conjugate_pr
%        'qf'  - ...conjugate_fr
%        'n' - minimization without derivatives, gsl_multimin_fminimizer_...
%        'nn'  - ...nmsimplex (default)
%        'n2'  - ...nmsimplex2
%        'nr'  - ...nmsimplex2rand
%  Initial parameters of optimization methods:
%    OPT.STEP    - 'step_size' for fdfminimizer_set, fminimizer_set (scalar) 
%    OPT.TOL     - 'tol' for fdfminimizer_set, fminimizer_set, 
%  Stopping parameters:  
%    OPT.EPSREL, OPT.EPSABS - 'gsl_multifit_test_delta' stopping criterion
%    OPT.EPSGRAD - 'gsl_multi..._test_gradient' stopping criterion
%  Advanced parameters:
%    OPT.REGGAMMA - regularization parameter for gamma, absolute
%
% PHI - transformation matrix (default I)
%
% Output arguments:
%
% XH   - low rank certificate: S(Ph) * PHI * [XH; -I] = 0
% INFO - convergence information: 
%   INFO.ITER - number of iterations
%   INFO.TIME - execution time
%   INFO.FMIN = (||P - Ph||_2)^2
% VXH  - asymptotic covariance matrix of vec(XH)
% PH   - Ph vector (approximation)
%
% Note: it is required that length(P) > size(S(P), 1) * (size(S(P), 2) - R)
