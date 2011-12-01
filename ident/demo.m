% DEMO - Demo file for identification by structured total least squares.

clear all
close all
echo off

% --- Zero lag corresponds to linear static model ---

l = 0;  % lag 
m = 2;  % # inputs
p = 3;  % # outputs
T = 20; % length of the data sequence

% Generate data
n    = l * p;           % order
sys0 = ss([], [], [], rand(p, m), 1);     % true system
u0   = randn(T, m);    % true input 
y0   = lsim(sys0, u0, 1:T); % true output
w0   = [u0 y0];       % true data  
w    = w0; % the given data is exact

% Identify the system
[sysh, info, wh, xini] = ident(w, m, l);
info

% Verify the results
err_sys = norm(sys0 - sysh)
err_w   = norm(w0 - wh)

pause

% --- Exact system identification ---

l = 2;    % lag
m = 2;    % # inputs
p = 3;    % # outputs
T = 100;  % length of the data sequence

% Generate data
n = l * p;             % order
sys0 = drss_(n, p, m); % true system
u0 = randn(T, m);     % true input 
y0 = lsim(sys0, u0, 1:T); % true output
w0 = [u0 y0];        % true data  
w  = w0;   % the given data is exact

% Identify the system
[sysh,info,wh,xini] = ident(w,m,l);
info

% Verify the results
err_sys = norm(sys0 - sysh)
err_w   = norm(w0 - wh, 'fro')

pause

% --- Errors-in-variables identification --- 

% Perturb the "true" data w0 with noise
w = w0 + 0.5 * randn(T, m + p);

% Identify the system
[sysh, info, wh, xini] = ident(w, m, l);
info

% Verify the results
err_data = norm(w0 - w , 'fro') / norm(w0, 'fro')
err_appr = norm(w0 - wh, 'fro') / norm(w0, 'fro')

% Plot "true", "noisy", and approximating input and output
echo off
t = 1:20;
for i = 1:m + p
  figure
  plot(t, w0(t, i), '-r', 'linewidth', 2), hold on
  plot(t, w(t, i) , ':k', 'linewidth', 2)
  plot(t, wh(t, i), '--b', 'linewidth', 2)
  set(gca, 'fontsize', 15)
  xlabel('time'), ylabel(sprintf('w_%d', i)), title('EIV identification')
  legend('true', 'data', 'appr')
end
echo on

pause

% --- Output error identification ---

% Perturb the "true" output y0 with noise
w    = w0 + 0.75 * [zeros(T, m) randn(T, p)];

% Identify the system
opt.exct = 1:m; % take into account that the input is exact
[sysh, info, wh, xini] = ident(w, m, l, opt);
info

% Verify the results
error_data = norm(w0 - w , 'fro') / norm(w0, 'fro')
error_appr = norm(w0 - wh, 'fro') / norm(w0, 'fro')

% Plot "true", "noisy", and approximating output 
echo off
t = 1:20;
for i = 1:p
  figure
  plot(t, w0(t, m + i), '-r', 'linewidth', 2), hold on
  plot(t, w(t, m + i) , ':k', 'linewidth', 2)
  plot(t, wh(t, m + i), '--b', 'linewidth', 2)
  set(gca, 'fontsize', 15)
  xlabel('time'), ylabel(sprintf('y_%d', i)), title('OE identification')
  legend('true', 'data', 'appr')
end
echo on

pause

% --- #inputs = 0 corrsponds to autonomous system identification ---

% Generate data
T     = 20; % Length of the data sequence
xini0 = randn(n, 1);
y0    = initial(sys0, xini0, T); y0 = y0(1:T, :);
w0    = y0;
% Perturb the "true" data w0 with noise
w     = w0 + 0.25 * randn(T, p); % data for the identification

% Approximate the given data by a system with complexity (0, l)
[sysh, info, wh, xinih] = ident(w, 0, l);
info

% Verify the results
error_data = norm(w0 - w , 'fro') / norm(w0, 'fro')
error_appr = norm(w0 - wh, 'fro') / norm(w0, 'fro')

% Plot "true", "noisy", and approximating sequences
echo off
t = 1:20;
for i = 1:p
  figure
  plot(t, w0(t, i), '-r', 'linewidth', 2), hold on
  plot(t, w(t, i) , ':k', 'linewidth', 2)
  plot(t, wh(t, i), '--b', 'linewidth', 2)
  set(gca, 'fontsize', 15)
  xlabel('time'), ylabel(sprintf('y_%d', i)), title('Autonomous system identification')
  legend('true', 'data', 'appr.')
end
echo on

pause

% --- Identification from step response observations ---

% Generate data
T   = 150;  % length of the data sequence
y0   = step(sys0, T); y0 = y0(1:T, :, :);
% Perturb the "true" data y0 with noise
y    = y0 + 0.25 * randn(T, p, m);

% Construct the inputs
u0   = zeros(T, m, m);
for i = 1:m
  u0(:, i, i) = ones(T, 1);
end
% The input/output data is
w = [u0 y];
% Proceed it with l zeros for the zero initial conditions
wext = [zeros(l, m + p, m); w];

% Identify the system from the extended I/0 data
opt.exct = [1:m]; % take into account that the input is exact
[sysh, info, whext, xini] = ident(wext, m, l, opt);
info

% Remove the trailing part, setting the initial conditions
wh = whext(l + 1:end, :, :);

% Verify the results
w0 = [u0 y0];
error_data = norm(w0(:) - w(:) ) / norm(w0(:))
error_appr = norm(w0(:) - wh(:)) / norm(w0(:))

% Plot "true", given, and approximating data sequences 
echo off
t = 1:20;
for i = 1:m
  for j = 1:p
    figure
    plot(t, w0(t, m + j, i), '-r', 'linewidth', 2), hold on
    plot(t, w(t, m + j, i) , ':k', 'linewidth', 2)
    plot(t, wh(t, m + j, i), '--b', 'linewidth', 2)
    set(gca, 'fontsize', 15)
    title('Identification from step response')
    xlabel('time'), ylabel(sprintf('s_{%d, %d}', i, j))
    legend('true', 'data', 'appr')
  end
end
echo on

pause

% --- Identification from impulse response observations ---
% ---       ( a.k.a. approximate realization )          ---

% Generate data
T  = 50; % length of the data sequence
y0 = impulse(sys0, T); y0 = y0(1:T, :, :);
% Perturb the "true" data y0 with noise
y  = y0 + 0.25 * randn(T, p, m);

% Version 1: solution by input/output identification.

% Construct the inputs
u0   = zeros(T, m, m);  u0(1, :, :) = eye(m);

% The given input/output data is
w = [u0 y];

% Proceed it with l zeros for the zero initial conditions
wext = [zeros(l, m + p, m); w];

% Identify the system from the extended I/0 data
opt.exct = [1:m]; % take into account that the input is exact
[sysh1, info1, whext, xini1] = ident(wext, m, l, opt);
info1 = info1

% Remove the trailing part, setting the initial conditions
wh1 = whext(l + 1:end, :, :);

% Version 2: solution by output-only identification.

% Identify an autonomous system
[sysh2, info2, yh2, xini2] = ident(y(2:end, :, :), 0, l);
yh2 = [y(1, :, :); yh2];
info2 = info2

% Recover the I/O system
sysh2 = ss(sysh2.a, xini2, sysh2.c, reshape(yh2(1, :, :), p, m), -1);
wh2 = [u0 yh2];

% Verify the results

w0 = [u0 y0];
error_data  = norm(w0(:)-w(:))/norm(w0(:))
error_appr1 = norm(w0(:)-wh1(:))/norm(w0(:))
error_appr2 = norm(w0(:)-wh2(:))/norm(w0(:))

% Plot "true", given, and approximating data sequences 
echo off
t = 1:20;
for i = 1:m
  for j = 1:p
    figure
    plot(t, w0(t, m + j, i), '-r', 'linewidth', 2), hold on
    plot(t, w(t, m + j, i) , ':k', 'linewidth', 2)
    plot(t, wh1(t, m + j, i), '--b', 'linewidth', 2)
    plot(t, wh2(t, m + j, i), '-.b', 'linewidth', 2)
    set(gca, 'fontsize', 15)
    title('Identification from impulse response')
    xlabel('time')
    ylabel(sprintf('h_{%d, %d}', i, j))
    legend('true', 'data', 'appr. 1', 'appr. 2')
  end
end
echo on

pause

% --- Finite-time l2 model reduction ---

l  = 10; % lag of the original (high order) system
lr = 1;  % lag of the reduced system
m  = 2;  % # inputs
p  = 2;  % # outputs

% High order system
sys = drss_(p*l, p, m);

% Simulate impulse response
h = impulse(sys); % determine automatically T
T = size(h, 1);    % time horizon T

% Find reduced model
[sysr, info, hr, xini] = ident(h(2:end, :, :), 0, lr);
hr = [h(1, :, :); hr];
sysr = ss(sysr.a, xini, sysr.c, reshape(hr(1, :, :), p, m), -1);
info

% Relative Hinf and H2 norms of the error system
h2_err  = norm(sys-sysr, 2) / norm(sys, 2)
hi_err  = norm(sys-sysr, 'inf') / norm(sys, 'inf')

% Plot the impulse responses of the original and appr. systems
echo off
t = 1:T;
for i = 1:m
  for j = 1:p
    figure
    plot(t, h(t, j, i), '-r', 'linewidth', 2), hold on
    plot(t, hr(t, j, i), '--b', 'linewidth', 2), hold on
    set(gca, 'fontsize', 15)
    title('Finite-time l2 model reduction')
    xlabel('time')
    ylabel(sprintf('h_{%d, %d}', i, j))
    legend('given', 'appr.', 'Location', 'Best')
  end
end
echo on

pause

% --- The misfit is generically independent of the input/output partitioning ---

l = 1;  % lag 
m = 2;  % # inputs
p = 1;  % # outputs
T = 20; % time horizon

% Generate data
sys0 = drss(p*l, p, m);
u0   = randn(T, m);
y0   = lsim(sys0, u0, 1:T);
w0   = [u0 y0];

% Identify the system from the original exact data
w1   = w0;
[sysh1, info1, wh1, xini1] = ident(w1, m, l);
info1

% Identify the system from exact data with randomly permuted variables
perm = randperm(m + p)
w2   = w1(:, perm);
[sysh2, info2, wh2, xini2] = ident(w2, m, l);
info2

% Compare the results
norm(wh1(:, perm) - wh2)

% --- End of demo --- 