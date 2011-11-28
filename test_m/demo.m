% DEMO - Demo file for structured total least squares

opt.maxiter = 100;
opt.epsrel  = 1e-7; 
opt.epsabs  = 0; 
opt.epsgrad = 1e-5; 
opt.disp    = 'disp';
opt.maxiter = 400;

addpath support
addpath stln
%addpath matlab_version

clear all
rand('state',0)
randn('state',0)
format long
echo on

% ===================================================
% ===================================================
% ======                                        =====
% ======              STLS demo                 =====
% ======                                        =====
% ===================================================
% ===================================================
%
%
% ===================================================
% ====== First we solve a least squares problem =====
% ===================================================
%
% Define dimensions and generate random data
m = 100; n = 5; d = 1;
a = rand(m, n);
b = rand(m, d);

% Find the LS estimate by Matlab's \
tic, x_ls = a \ b; t_ls = toc
% x_ls(1,1:d)

% Define and solve the LS problem as a (very special) STLS problem
s_ls.a = [4 n 1; 3 d 1]; s_ls.m = m; s_ls.d = d; s_ls.k = 1;
pause
at = a'; bt = b';
%tic, 
[x_stls, i_stls] = stls([], [], s_ls, [], [], [at(:); bt(:)]); 
%t_stls = toc
% x_stls(1,1:d)

% Press any key to continue 
pause

% ===================================================
% ======      Next we solve a TLS problem       =====
% ===================================================
%
% The data is a,b used above.
%
% Solve the TLS problem via SVD

a0 = rand(m, n); a = a0 + 0.2 * randn(m, n);
x0 = rand(n, d);
b0 = a0 * x0;  b = b0 + 0.2 * randn(m, d);

tic, x_tls = tls(a,b); t_tls = toc
x_tls(1,1:d)

% Define and solve the TLS problem as an STLS problem
s_tls = [3 n + d 1];
%tic, 
[x_stls, i_stls] = stls(a, b, s_tls); 
%t_stls = toc
x_stls(1, 1:d)

% Press any key to continue
pause

%% ============================================================
%% ======      Next we solve a mixed LS-TLS problem       =====
%% ============================================================
%%
%% The data is a,b used above.
%%
%n1 = 1; % # of column of a1, where a =: [a1 a2] with a1 exact

%% Solve the LS-TLS problem via exact algorithm
%tic, x_lstls = ls_tls(a(:,1:n1),a(:,n1+1:end),b); t_lstls = toc
%x_lstls(1,1:d)

%% Define and solve the LS-TLS problem as an STLS problem
%s_lstls = [4 n1 1; 3 n+d-n1 1];
%tic, [x_stls,i_stls] = stls(a,b,s_lstls); t_stls = toc
%x_stls(1,1:d)

% Press any key to continue
%pause

% ===================================================
% ======  Next we solve a deconvolution problem =====
% ===================================================
%
% b0 = conv(p_a0,x0), p_a = p_a0 + noise, b = b0 + noise
% Problem: given a and b, estimate x0

m = 200; % length(p_a0)
n = 2;   % length(b0)

% Generate true data: p_a0 and b0
p_a0 = rand(n + m - 1);
a0   = toeplitz(p_a0(n:n+m-1),p_a0(n:-1:1));
x0   = rand(n,1);
b0   = a0*x0;

% Add noise: p_a = p_a0 + noise, b = b0 + noise
v_n  = 0.25; % noise level
p_a  = p_a0 + v_n * randn(n+m-1);
a    = toeplitz(p_a(n:n+m-1),p_a(n:-1:1));
b    = b0 + v_n * randn(m,1);

% Ignore the structure and estimate via LS and TLS
tic, xh_ls  = a\b; t_ls = toc
tic, xh_tls = tls(a,b); t_ls = toc

% Define the structure and solve the deconvolution problem via STLS
s = [1 n 1; 3 1 1];
tic, [xh_stls,i_stls] = stls(a,b,s); t_stls = toc
xh_stls(1:2)'
i_stls.fmin % value of the cost function at xh_stls

% Solve via an alternative STLS method
tic, xh_stln = faststln1(a,b, xh_ls); t_stln = toc
xh_stln(1:2)'
cost1(xh_stln,a,b,s) % value of the cost function at xh_stln

% Compare the relative errors of estimation 
e_ls   = norm(xh_ls-x0)/norm(x0); disp(e_ls)
e_tls  = norm(xh_tls-x0)/norm(x0); disp(e_tls)
e_stls = norm(xh_stls-x0)/norm(x0); disp(e_stls)

% Press any key to continue
pause

% ==================================================================
% ====== Next we solve a Hankel low rank approximation problem =====
% ==================================================================
%
% Generate data
%np = 12;                      % number of parameters
%p0 = (1:np)';                 % true value of the parameter vector
%p  = p0 + [5; zeros(np-1,1)]; % add disturbance
%c = hankel(p(1:10),p(10:np));
%a = c(:,1:2); b = c(:,3);

% Define the structure and solve the problem via STLS 
%s  = [2 3 1];
%tic, [xh_stls,i_stls] = stls(a,b,s); t_stls = toc
%[dp, ch_stls ] = corr( xh_stls, [a b], s );
%i_stls.fmin % value of the cost function at xh_stls

% Check the result
%svd_ch_stls = svd(ch_stls); disp(svd_ch_stls)
%norm_res = norm(ch_stls * [xh_stls; -1]); disp(norm_res)

% Solve via an alternative STLS method
%c = [a b]; c = fliplr(c);
%tic, xh_stln = faststln2(c(:,1:2),c(:,3)); t_stln = toc
%x_ext = [xh_stln; -1]; x_ext = flipud(x_ext); xh_stln = -x_ext(1:2)/x_ext(3);
%xh_stln(1:2)';
%cost1(xh_stln,a,b,s) % value of the cost function at xh_stln

% Press any key to continue
%pause

% % =====================================================================
% % ====== Next we solve a transfer function identification problem =====
% % =====================================================================
% 
% % True model
% n   = 3;
% num = 0.151*[1 0.9 0.49 0.145];
% den = [1 -1.2 0.81 -0.27];
% 
% % True data
% t  = 50;
% u0 = randn(t,1); 
% [y0,x0] = dlsim(num,den,u0);
% 
% % Noisy data
% v_noise = .01;
% y = y0 + v_noise * randn(t,1);
% u = u0 + v_noise * randn(t,1);
% 
% % Press any key to continue
% pause
% 
% % Define a and b
% m = length(y)-n;
% a = [ hankel(u(1:m),u(m:end)) hankel(y(1:m),y(m:end)) ];
% b = a(:,end);
% a(:,end) = [];
% 
% % Ignore the structure and solve the identification problem via LS and TLS
% tic, xh_ls  = a\b; t_ls = toc
% tic, xh_tls = tls(a,b); t_tls = toc
% 
% % Define the structure and solve the identification problem via STLS
% s = [2 n+1 1; 2 n+1 1];
% tic, [xh_stls,i_stls] = stls(a,b,s); t_stls = toc
% 
% % Extract the estimates 
% num_ls = fliplr(xh_ls(1:n+1)'); den_ls = [1 fliplr(-xh_ls(n+2:end)')];
% num_tls = fliplr(xh_tls(1:n+1)'); den_tls = [1 fliplr(-xh_tls(n+2:end)')];
% num_stls = fliplr(xh_stls(1:n+1)'); den_stls = [1 fliplr(-xh_stls(n+2:end)')];
% 
% % Compare the relative errors of estimation
% e_ls = norm([num - num_ls, den - den_ls]) / norm([num,den]); disp(e_ls)
% e_tls = norm([num - num_tls, den - den_tls]) / norm([num,den]); disp(e_tls)
% e_stls = norm([num - num_stls, den - den_stls]) / norm([num,den]); disp(e_stls)

% ===================================================
% ======              End of demo              ======
% ===================================================
echo off
