function [M, wh, xini] = misfit(w, sys, exct)
% MISFIT - Global Total Least Squares misfit ||W - Wh||_F. 
% 
% [M, wh, xini] = misfit(w, sys, exct)
% 
% Inputs:
% W   - given a set of time series, stored in an arry with dimensions
%                 #samples x #variables x #time series            (*)
%       in case of equal number of samples or a cell array with entries 
%       of the type (*) in case of different number of samples
% SYS  - state-space system with M inputs, P := size(W, 2) - M outputs, 
%        and order N := L * P, where L is an integer
% EXCT - vector of indices for exact variables (default [])
%
% Outputs:
% M    - misfit ||W - WH||_F
% WH   - optimal approximating time series
% XINI - initial condition, under which WH is obtained

%% Constants
[p, m] = size(sys);
n      = size(sys, 'order');
l      = n / p;

%% Call ident with initial approximation sys and zero iterations
opt.sys0 = sys; opt.maxiter = 0; 
if nargin == 3, opt.exct = exct; end
[sysh, info, wh, xini] = ident(w, m, l, opt); M = info.M;