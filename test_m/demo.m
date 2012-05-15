% DEMO - Demo file for structured low rank approximation
clear all, rand('state', 0), randn('state', 0), addpath '..';

opt.disp = 'iter';

%% Least squares problem

% Define dimensions and generate random data
m = 100; n = 5; d = 2; sigma = 0.1;
a0 = rand(m, n); x0 = rand(n, d); b0 = a0 * x0;
a = a0; b = b0 + sigma * randn(m, d);

% Find the LS estimate by Matlab's \
tic, x_ls = (a \ b)'; t_ls = toc

% Define and solve the LS problem as a (very special) SLRA problem
s_ls.m = ones(1, n+d);
s_ls.n = m; 
s_ls.w = (diag([Inf Inf Inf Inf Inf 1 1]) * ones(n+d, m))';
p = [a b];
tic, [p_slra, i_slra] = slra(p, s_ls, n, opt); t_slra = toc
error = i_slra.Rh(:,1:end-d) - x_ls

pause
%% Low rank approximation problem

% Define dimensions and generate random data
%m = 100; n = 5; d = 2; sigma = 0.1;
D0 = [a0 b0]; %rand(m, n) * rand(n, n + d); 
D  = D0 + sigma * randn(m, n + d);

tic, x_lra = tls(D(:, 1:n), D(:, n + 1:end))'; t_lra = toc, 

% Define and solve the LRA problem as an SLRA problem
s_lra.m = ones(1, n+d);
s_lra.n = m; p = D(:);
opt.Rini = i_slra.Rh;
tic, [p_slra, i_slra] = slra(p, s_lra, n, opt); t_slra = toc
error = i_slra.Rh(:,1:end-d) - x_lra

pause
%% Deconvolution problem

m = 200; % length(p_a0)
n = 2;   % length(b0)
sigma = 0.25; % noise level

% Generate true data: p_a0 and b0
p_a0 = rand(n + m - 1, 1);
a0   = toeplitz(p_a0(n:n + m - 1), p_a0(n:-1:1));
x0   = rand(n, 1);
b0   = a0 * x0;

% Add noise: p_a = p_a0 + noise, b = b0 + noise
p_a = p_a0 + sigma * randn(n + m - 1, 1);
a   = toeplitz(p_a(n:n + m - 1), p_a(n:-1:1));
b   = b0 + sigma * randn(m, 1);

% Ignore the structure and estimate via LS and TLS
tic, xh_ls  = a \ b; t_ls = toc
tic, xh_tls = tls(a, b); t_tls = toc

% Solve the deconvolution problem via SLRA
s.m = [n 1];
s.n = 200;
s.phi = [0 1 0; 1 0 0; 0 0 1];
opt =  rmfield(opt, 'Rini');
tic, [p_slra, i_slra] = slra([p_a; b], s, n, opt); t_slra = toc

% Compare the relative errors of estimation 
e_ls   = norm(xh_ls - x0) / norm(x0)
e_tls  = norm(xh_tls - x0) / norm(x0)
e_slra = norm((i_slra.Rh(n:-1:1))' - x0) / norm(x0)

pause
%% Hankel structured low rank approximation problem

% Generate data
np = 50;                            % number of parameters
p0 = (1:np)'; n = 2                 % true value of the parameter vector
sigma = 0.1;                        % noise standard deviation
p  = p0 + sigma * randn(np, 1);          % add disturbance
c = hankel(p(1:np - n), p(np - n:np));
a = c(:, 1:n); b = c(:, n + 1);

% Define the structure and solve the problem via SLRA 
s.m = n+1;
s.n = np-n;
s.phi = eye(n+1);
tic, [p_slra, i_slra] = slra(p, s, n, opt); t_slra = toc
error_data = norm(p - p0)
error = norm(p_slra - p0)

% Make the same with a permutation matrix
I = eye(n + 1); perm = randperm(n + 1); 
s.phi = I(perm, :);
tic, [p_slra2, i_slra] = slra(p, s, n, opt); t_slra = toc
error_data = norm(p - p0)
error = norm(p_slra2 - p0)

norm(p_slra - p_slra2)

pause

%% Sylvester low-rank approximation 
% Example 1
clear all;
opt.disp = 'iter';
a = conv([4 2 1], [5 2]) + [0.05 0.03 0.04 0]
b = conv([4 2 1], [5 1]) + [0.04 0.02 0.01 0]
s.m = [2 2];
s.n = 5;
s.w = [Inf ones(1,length(a)) Inf Inf ones(1,length(b)) Inf];
[ph, info] = slra([0 a 0 0 b 0], s, 3, opt)


% Example 2 (from book)
opt.disp = 'iter';
a2 = conv([-1 1], [-5 1])
b2 = conv([-1.1 1], [-5.2 1])
s.m = [2 2];
s.n = 4;
s.w = [Inf ones(1,length(a2)) Inf Inf ones(1,length(b2)) Inf];
[ph, info] = slra([0 a2 0 0 b2 0], s, 3, opt)

pause
% Example 3 (3 polynomials)
opt.disp = 'iter';
a3 = conv([-1 1], [-5.1 1])
b3 = conv([-1.1 1], [-5 1])
c3 = conv([-1 1], [-5.2 1])
p = [0 a3 0 b3 0 zeros(1,5) c3 zeros(1,5)]
s = struct('m', [2 6], 'n', 8);
s.w = [Inf ones(1,length(a3)) Inf ones(1,length(b3)) Inf Inf * ones(1, 5) ones(1,length(c3)) Inf * ones(1, 5)];
s.phi = [eye(2) zeros(2,6); zeros(2,2) eye(2) zeros(2,4); zeros(2,6) eye(2)];
%opt.reggamma = 1;
%opt.reggamma = 10000;
[ph, info] = slra(p, s, 5, opt)





