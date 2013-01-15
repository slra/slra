% MISFIT - Global Total Least Squares misfit ||W - Wh||_F. 
% 
% [M, wh, xini] = misfit(w, sys, opt)
% 
% Inputs:
% W   - given a set of time series, stored in an arry with dimensions
%                 #samples x #variables x #time series            (*)
%       in case of equal number of samples or a cell array with entries 
%       of the type (*) in case of different number of samples
% SYS - state-space system with M inputs, P := size(W, 2) - M outputs, 
%        and order N := L * P, where L is an integer
% OPT - options 
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.WINI = 0 specifies zero initial conditions (default [])
%
% Outputs:
% M    - misfit ||W - WH||_F
% WH   - optimal approximating time series
% XINI - initial condition, under which WH is obtained
function [M, wh, xini] = misfit(w, sysh, varargin)
[p, m] = size(sysh); n = size(sysh, 'order'); ell = n / p;
opt.sys0 = sysh; opt.maxiter = 0; opt.disp = 'off';
[sysh, info, wh, xini] = ident(w, m, ell, varargin{:}, opt); M = info.fmin;
