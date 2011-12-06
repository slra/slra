% SLRA - Structured Low Rank Approximation.
%
% Solves the structured system of equations A*X = B with size(A,1) > size(A,2).
% The augmented data matrix C = [A B] is of the form C = [C1 ... Cq], where 
% the blocks Ci are (block) Toeplitz or Hankel structured, unstructured, or 
% constant. C depends affinely on a vector of parameters P, i.e., C = C(P). 
% The STLS problem is: min_{Xh,Ph} (P-Ph)'*(P-Ph) s.t. C(Ph)*[Xh;-I] = 0.
%
% [ xh, info, vxh, fh ] = stls( p, s, r, x0, opt )
%
% Input arguments:
%
% P - parameter vector.
% S - structure specification, a matrix with q rows and 1-4 columns:
%   S(i,1) - number of subblocks in each block row of the i-th block (put 1 if you need unstructured matrix)
%   S(i,2) - number of columns in each subblock (by default is 1)
%   S(i,3) - = 1 if the block is exact,  = 0 otherwise (default)
%   S(i,4) - = 1 if the block is Toeplitz,  = 0 if it is Hankel (default)
% For block-Toeplitz/Hankel structured problem with K x nu(i) size blocks, 
% S is a structure with fiels K and A, where S.A is the matrix described 
% above and S.K = K.
% R - rank of the approximation (size(A,2), by default it is equal to sum(S(:,2))) -1 
% X0   - initial approximation (default TLS).
% OPT  - optimization parameters, structure with fields: 
%   OPT.MAXITER - maximum number of iterations, 
%   OPT.EPSREL, OPT.EPSABS, OPT.EPSGRAD - convergence tolerances, and 
%   OPT.REGGAMMA - regularization parameter for gamma, absolute (to be changed)
%   OPT.DISP - level of display ['notify','final','iter',off] (default 'notify').
% Exit condition: #iterations >= MAXITER, |xi(k+1)-xi(k)| < EPSABS + EPSREL*|xi(k+1)|
% for all i, where x(k) is the approximation on the i-th step, or ||f'|| < EPSGRAD, 
% where f' is the gradient of the cost function.
%
% Output arguments:
%
% XH   - STLS estimate.
% INFO - information on exit, structure with fields ITER, TIME, and FMIN:
%   INFO.ITER - number of iterations for convergence,
%   INFO.TIME - execution time,
%   INFO.FMIN - value of the cost function.
% VXH  - assymptotic covariance matrix of the estimate XH.
%
% Note: The program can not treat the case length(P) < size(A,1) * size(B,2).

% Author: Ivan Markovsky, Last modified: November 2004.
%
% Reference: I. Markovsky and S. Van Huffel "High-performance numerical algorithms 
% and software for structured total least squares", Journal of Computational and 
% Applied Mathematics, 2005
