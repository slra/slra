%% REG_SLRA - structured low rank approximation using factorization approach
%
%  Approximates a given structured (m x n) matrix S(p)
%  as a product PL of two factors P(m x r), L:(r x n),
%  such that PL is again a structured matrix.
%
%% Syntax
%  [ph, info] = slra_reg(p, s, r, opt)
%
%% Input 
% p       - structure parameter vector
% s       - structure specification s
% - s.tts - integer matrix, see the documentation of slra function
% - s.S0  - S0 [optional; default: zero matrix]
% r       - rank of the approximation
% opt     - options [optional]
% - opt.w - vector of weights [default: vector of ones]
% - opt.P_init - initial approximation for P [default: based on SVD]
% - opt.lambda_init - initial value of lambda [default: 1]
% - opt.lambda_max - maximal value of lambda [default: 1e+14]
% - opt.max_inner_iter - maximal # inner steps [default: 20];
%
%% Output
% ph        - approximation of p, corresponding to low-rank matrix
% info      - other outputs
% - info.P  - P matrix
% - info.L  - L matrix
% - info.lambda - last value of lambda
% - info.outer_iterations - # outer iterations
% - info.inner_iterations - # inner iterations for each outer iteration
% - info.fmin = Mslra_ext(info.Rh, s.tts, p, opt.w);
% - info.error_constraints = norm(PL_1 - PL_2, 'fro')^2 / norm(PL_1, 'fro')^2;
%
%% Reference
% M. Ishteva, K. Usevich, and I. Markovsky.
% Regularized structured low-rank approximation with applications. 
% SIAM J. Matrix Anal. Appl., 35(3):1180-1204, 2014.
%
%% See also
%   slra
function [ph, info] = reg_slra(p, s, r, opt)

% set parameters
[s.m, s.n] = size(s.tts);
s.S = zeros(s.m * s.n, length(p)); % construct S
for i = 1:length(p)
    s.S(s.tts==i,i) = 1;
end
if ~isfield(s, 'S0') || isempty(s.S0) % construct S0
    s.S0 = zeros(s.m * s.n, 1);
end
if ~exist('opt', 'var') % check options
    opt = [];
end
if ~isfield(opt, 'w') || isempty(opt.w) % weights
    opt.w = ones(length(p),1);
    M = diag(opt.w);
else
    M = sqrtm(diag(opt.w));
end
if ~isfield(opt, 'P_init') || isempty(opt.P_init) % initial approximation
    [opt.P_init, ~, ~] = svds(reshape(s.S0 + s.S * p, s.m, s.n), r);
end
if ~isfield(opt, 'lambda_init') || isempty(opt.lambda_init) % initial lambda
    lambda_init = 1;
end
if ~isfield(opt, 'lambda_max') || isempty(opt.lambda_max) % maximal lambda
    lambda_max = 1e+14;
end
if ~isfield(opt, 'max_inner_iter') || isempty(opt.max_inner_iter) % maximal # inner steps
    max_inner_iter = 20;
end

lambda_small_step = 1.5;
lambda_large_step = 10;
lambda = lambda_init;
P = opt.P_init;
P_old = randn(size(P));
L = zeros(r, s.n);

Pi_S = (s.S' * s.S) \ s.S';
rhs = sparse([M * p;  sqrt(lambda) * s.S0(:)]);
lf = sparse([M * Pi_S; sqrt(lambda) * (eye(s.m * s.n) - s.S * Pi_S)]);

% main loop
i = zeros(100, 1);
j = 0;
while lambda < lambda_max
    j = j + 1;
    % inner loop
    k = 0;
    while subspace(P, P_old) >  1e-3 / sqrt(lambda) && k < max_inner_iter
        k = k + 1;
        P_old = P;
        L = full(reshape((lf * sparse(kron(eye(s.n), P))) \ rhs, r, s.n));
        P = full(reshape((lf * sparse(kron(L', eye(s.m)))) \ rhs, s.m, r));
    end 
    i(j) = k;
    P_old = randn(size(P));
    % update lambda
    if k < 10
        lambda = lambda * lambda_large_step;
    else
        lambda = lambda * lambda_small_step;
    end
    rhs = sparse([M * p;  sqrt(lambda) * s.S0(:)]);
    lf = sparse([M * Pi_S; sqrt(lambda) * (eye(s.m * s.n) - s.S * Pi_S)]);
end

% output
info.P = P;
info.L = L;
info.lambda = lambda;
info.outer_iterations = j;
info.inner_iterations = i(1:j);
info.fmin = Mslra_ext(null(full(P)')', s.tts, p, opt.w);
PL_1 = P * L;
ph = Pi_S * PL_1(:);
PL_2 = reshape(s.S0 + s.S  * ph, s.m, s.n);
info.error_constraints = norm(PL_1 - PL_2, 'fro')^2 / norm(PL_1, 'fro')^2;