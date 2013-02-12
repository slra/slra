% IDENT_MISFIT - orthogonal distance from W to an LTI system SYS
% M = minimum over Wh  ||W - Wh||_F^2 s.t. Wh is a trajectory of SYS
% [M, wh, xini] = ident_misfit(w, sys, opt)
% 
% W   - time series or a set of time series (see "help ident")
% SYS - state-space system with M inputs, P := size(W, 2) - M outputs
%        and order N := L * P, where L is an integer
% OPT - options 
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.WINI = 0 specifies zero initial conditions (default [])
% M    - misfit ||W - WH||_F^2
% WH   - optimal approximating time series
% XINI - initial condition under which WH is obtained by SYSH
function [M, wh, xini] = ident_misfit(w, sysh, opt)
[p, m] = size(sysh); n = size(sysh, 'order'); ell = n / p;
opt.sys0 = sysh; opt.maxiter = 0; opt.disp = 'off';
[sysh, info, wh, xini] = ident(w, m, ell, opt); M = info.M;
