clear all
%addpath ~/mfiles/ident
%addpath ~/springer

l = 2;    % lag
m = 8;    % # inputs
p = 12;    % # outputs
n = l * p;             % order
sys0 = drss(n, p, m); % true system

% Generate data
T  = 50; % length of the data sequence
y1 = shiftdim(impulse(sys0, T + 1),1);
y0 = y1(:,:,2:T+1)

plot(reshape(y0(1,1,:),50,1));
%for i = 1:m
%    y0i = impulse(sys0, i, T + 1); y0i = y0i(:, 2:(T + 1));
%    y0(:, i, :) = reshape(y0i(:), p, 1, T);
%end
% Perturb the "true" data y0 with noise
y  = y0 + 0.01 * randn(p,m,T);


h0 = blkhankel(y0, 3);

svd(h0);
svd(h);


h = blkhankel(y, 3);
h = h'
a = h(:,1:n);
b = h(:,(n+1):(3*p));


h0 = h0'
a0 = h0(:,1:n);
b0 = h0(:,(n+1):(3*p));

plot([reshape(y0(1,1,:),50,1) reshape(y(1,1,:),50,1)]);

