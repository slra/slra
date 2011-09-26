clear all
addpath ~/mfiles/ident
addpath ~/springer

l = 2;    % lag
m = 8;    % # inputs
p = 12;    % # outputs
n = l * p;             % order
sys0 = drss_(n, p, m); % true system

% Generate data
T  = 50; % length of the data sequence
y0 = zeros(p, m, T);
for i = 1:m
    y0i = impulse(sys0, i, T + 1); y0i = y0i(:, 2:(T + 1));
    y0(:, i, :) = reshape(y0i(:), p, 1, T);
end
% Perturb the "true" data y0 with noise
y  = y0 + 0.25 * randn(p,m,T);

