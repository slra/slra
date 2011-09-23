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

% Version 1: solution by input/output identification.

% Construct the inputs
u0   = zeros(T,m,m);  u0(1,:,:) = eye(m);

% The given input/output data is
w = [u0 y];

% Proceed it with l zeros for the zero initial conditions
wext = [zeros(l,m+p,m); w];

% Identify the system from the extended I/0 data
opt.exct = [1:m]; % take into account that the input is exact
[sysh1,info1,whext,xini1] = stlsident(wext,m,l,opt);
info1 = info1

% Remove the trailing part, setting the initial conditions
wh1 = whext(l+1:end,:,:);

% Version 2: solution by output-only identification.

% Identify an autonomous system
[sysh2,info2,yh2,xini2] = stlsident(y(2:end,:,:),0,l);
yh2 = [y(1,:,:); yh2];
info2 = info2

% Recover the I/O system
sysh2 = ss(sysh2.a,xini2,sysh2.c,reshape(yh2(1,:,:),p,m),-1);
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
    plot(t,w0(t,m+j,i),'-r','linewidth',2), hold on
    plot(t,w(t,m+j,i) ,':k','linewidth',2)
    plot(t,wh1(t,m+j,i),'--b','linewidth',2)
    plot(t,wh2(t,m+j,i),'-.b','linewidth',2)
    set(gca,'fontsize',15)
    title('Identification from impulse response')
    xlabel('time')
    ylabel(sprintf('h_{%d,%d}',i,j))
    legend('true','data','appr. 1','appr. 2')
  end
end
echo on
