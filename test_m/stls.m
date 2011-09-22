% STLS - Structured Total Least Squares approximation.
%
% Solves the structured system of equations A*X = B with size(A,1) > size(A,2).
% The augmented data matrix C = [A B] is of the form C = [C1 ... Cq], where 
% the blocks Ci are (block) Toeplitz or Hankel structured, unstructured, or 
% constant. C depends affinely on a vector of parameters P, i.e., C = C(P). 
% The STLS problem is: min_{Xh,Ph} (P-Ph)'*(P-Ph) s.t. C(Ph)*[Xh;-I] = 0.
%
% [ xh, info, vxh ] = stls( a, b, s, x0, opt )
%
% Input arguments:
%
% A, B - data.
% S    - structure specification, a matrix with q rows and 3 columns:
%   S(i,1) - structure type [ 1 | 2 | 3 | 4 ] of the i-th block,
%     1 - Toeplitz, 2 - Hankel, 3 - Unstructured, 4 - Exact (noise free),
%   S(i,2) - number of columns of the i-th block,
%   S(i,3) - number of columns of the blocks in the i-th block.
% For block-Toeplitz/Hankel strucutred problem with K x nu(i) size blocks, 
% S is a structure with fiels K and A, where S.A is the q x 2 matrix described 
% above and S.K = K.
% X0   - initial approximation (default TLS).
% OPT  - optimization parameters, structure with fields: 
%   OPT.MAXITER - maximum number of iterations, 
%   OPT.EPSREL, OPT.EPSABS, OPT.EPSGRAD - convergence tolerances, and 
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
