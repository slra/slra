function [sysh, info, w, xini] = ident(w, m, l, opt)
% IDENT - Identification by structured low rank approximation
%
% [sysh, info, wh, xini] = ident(w, m, l, opt)
%
% Inputs:
% W   - given a set of time series, stored in an arry with dimensions
%                 #samples x #variables x #time series            (*)
%       in case of equal number of samples or a cell array with entries 
%       of the type (*) in case of different number of samples
% M   - input dimension
% L   - system lag
% OPT - options for the optimization algorithm:
%   OPT.EXCT - vector of indices for exact variables (default []);
%   OPT.SYS0 - initial approximation: an SS system with M inputs, 
%              P := size(W, 2) - M outputs, and order N : =L * P;
%   OPT.DISP    - level of display [off | iter | notify | final];
%   OPT.MAXITER - maximum number of iterations (default 100);
%   OPT.EPSREL, OPT.EPSABS, and OPT.EPSGRAD (default 1e-4)
%     - convergence tolerances (see the help of SLRA)
%
% Outputs:
% SYSH - input/state/output representation of the identified system
% INFO - information from the optimization solver:
%   INFO.M - misfit ||W - WH||_F,
%   INFO.TIME - execution time for the SLRA solver,
%   INFO.ITER - number of iterations
%   NOTE: INFO.ITER = OPT.MAXITER indicates lack of convergence
% W    - optimal approximating time series
% XINI - initial conditions under which WH are obtained

%% Constants
if ~iscell(w)
    T  = size(w, 1); % # of samples
    nw = size(w, 2); % # of variables
    N  = size(w, 3); % # of time series
else
    error('Cell array input is not implemented yet.')
end
p  = nw - m;     % # of outputs
n  = p * l;      % order of the system
siso = (m == 1) & (p == 1);
l1   = l + 1; 

%% Check for parameters out of range
if nargin < 3
    error('Not enough input arguments. (W, M, and L should be given.)')
end
if (~isreal(w))
    error('Complex data not supported.')
end
if (~isreal(m) | m < 0 | m > nw | (ceil(m) - m ~=  0))
    error('Invalid number of inputs M. (0 < M < size(W, 2) and integer.)')
end
if (~isreal(l) | l < 0 | (ceil(l) - l ~=  0))
    error('Invalid system lag L. (0 < L and integer.)')
end

%% Default options values
disp    = 'notify';
epsrel  = 1e-4;
epsabs  = 1e-4;
epsgrad = 1e-4;
maxiter = 10;
sys0    = [];  % to be chosen
exct    = [];  % no exact variables

%% User specified options
OPTNAMES = {'disp', 'epsrel', 'epsabs', 'epsgrad', ...
            'maxiter', 'sys0', 'exct'};
DISPNAMES = {'off', 'iter', 'notify', 'final'};
if (nargin > 3)
    if isstruct(opt)
        names = fieldnames(opt);
        I = length(names);
        for i = 1:I
            ind = strmatch(lower(names{i}), OPTNAMES, 'exact');
            if isempty(ind)
                warning(sprintf('OPT.%s is an invalid option. Ignored.', ...
                                upper(names{i})))
            else
                eval([OPTNAMES{ind} ' = opt.' names{i} ';']);
            end
        end
    else
        warning('OPT should be a structure. Ignored.')
    end
end
opt = [];
opt.epsrel  = epsrel;
opt.epsabs  = epsabs;
opt.epsgrad = epsgrad;
opt.maxiter = maxiter;
opt.disp    = disp;

%% Check for invalid option values
exct = unique(exct);
if any(1 > exct | nw < exct)
    error('An index for exact variable is out of range.')
end
if (length(exct) * T > m * T + n)
    error('Too many exact variables. Generically there is no solution.')
end

%% Detect trivial cases
if (m == nw) 
    warning('M = size(W, 2) => trivial solution SYSH = I, WH = W.');
    sysh      = ss([], [], [], eye(nw), 1);
    info.M    = 0;
    info.iter = 0;
    if nargout == 4
        xini = inistate(w, sysh);
    end
    return
end
if (T <= n) 
    warning('T < p * l => trivial solution.');
    a = diag(ones(n - 1, 1), -1);
    b = zeros(n, m);
    c = [w(:, m + 1:end)' zeros(p, T - n)];
    d = zeros(p, m);
    sysh      = ss(a, b, c, d, 1);
    info.M    = 0;
    info.iter = 0;
    if nargout == 4
        xini = inistate(w, sysh);
    end
    return
end

%% Initial approximation
if ~isempty(sys0)
    if ~isa(sys0, 'ss')
        warning('OPT.SYS0 not an SS object. Ignored.')
        sys0 = [];
    else
        [pp, mm] = size(sys0); 
        nn = size(sys0, 'order');
        if (mm ~= m) | (pp ~= p) | (nn ~= n)
            warning('OPT.SYS0 invalid. Ignored.')
            sys0 = [];
        end
    end
end
if isempty(sys0)
    %[a, b, c, d] = subid(w(:, m + 1:end), w(:, 1:m), l + 1, n, [], [], 1); 
    %sys0 = uy2ssbal(w(:, 1:m), w(:, m + 1:end), l, 5 * l); 
    %sys0 = ss(a, b, c, d, 1);
end

%% Structure specification
ne = length(exct);
nn = nw - ne;
noisy = 1:nw; noisy(exct) = []; % indices of noisy variables
if ne == 0
    struct = [2 nw * l1 nw]; par = vec(shiftdim(w, 1));
else % Define a block of exact variables
    struct = [4 ne * l1 ne; 2 nn * l1 nn]; 
    par = [vec(blkhankel(w(:, exct, :), l1));
           vec(shiftdim(w(:, noisy, :), 1))];
end
if N > 1
    struct.a = struct;
    struct.k = N;
end

%% Call the SLRA solver
x0 = sys2x(sys0, exct, noisy);
[x, info, cov, parh] = slra(par, struct, l1 * m + n, x0, opt);
info.M = sqrt(info.fmin); info = rmfield(info, 'fmin');

if nargout > 2
    if N == 1
        w(:, noisy) = reshape(parh(end - nn * T + 1:end), nn, T)';
    else % N > 1
        w(:, noisy, :) = shiftdim(reshape(parh(end - nn * T * N + 1:end), nn, N, T), 2);
    end
end

%% Find an I/S/O representation of the model
if m > 0 % treat separately I/O and output only cases
    R = [x', -eye(p)]; % kernel representation of SYSH
    if ne > 0
        R3e = reshape(R(:, 1:ne * l1), p, ne, l1);
        R3n = reshape(R(:, l1 * ne + 1:end), p, nn, l1);
        R3  = [R3e R3n];
        R3(:, [exct noisy], :) = R3; % restore the original ordering
        clear R3n R3e
    else % ne == 0
        R3 = reshape(R, p, nw, l1);
    end
    Q3 =  R3(:, 1:m, :); % extract P and Q
    P3 = -R3(:, m + 1:nw, :);
    a  = zeros(n); b  = zeros(n, m); c = [];
    if n > 0
        a(p + 1:end, 1:n - p) = eye(n - p);
        c  = [zeros(p, n - p) eye(p)]; 
    end
    d  = Q3(:, :, l1);
    ind_j = (n - p + 1):n;
    for i = 1:l
        ind_i = (i - 1) * p + 1:i * p;
        P3i = P3(:, :, i);
        a(ind_i, ind_j) = -P3i;
        b(ind_i, :) = Q3(:, :, i) - P3i * d;
    end
    sysh = ss(a, b, c, d, 1);
else % output only case 
    O    = null([x' -eye(p)]);
    ch   = O(1:p, :);
    ah   = O(1:end-p, :) \ O(p + 1:end, :);
    sysh = ss(ah, [], ch, [], 1);      
end

%% Find the initial conditions
if nargout > 3    
    if nargout > 3
        xini = inistate(w, sysh);
    else
        xini = zeros(n, N);
        for i = 1:N
            tmp  = w(1:l1, :, i)'; 
            xini(:, i) = O \ tmp(:);
        end
    end
end


%% VEC - vectorize a matrix
function a = vec(A), a = A(:);


%% BLKHANKEL - construct a block Hankel matrix
function H = blkhankel(w, i)

[T, nw, N]  = size(w); j = T - i + 1;
if j <= 0
  error('Not enough data.')
end
H = zeros(i * nw, j * N);
if (N > 1) % multiple time series => block matrix
  w = reshape(shiftdim(w, 1), nw, N * T);
  for ii = 1:i
    H((ii - 1) * nw + 1:ii * nw, :) = w(:, (ii - 1) * N + 1:(ii + j - 1) * N);
  end
else % simgle time series => block vector
  w = w';
  for ii = 1:i
    H((ii - 1) * nw + 1:ii * nw, :) = w(:, ii:ii + j - 1);
  end
end


%% SYS2X - Convert (A,B,C,D) to a parameter X for SLRA
function x = sys2x(sys, exct, noisy)

if ~isa(sys, 'ss')
    x = []; % use default for the slra solver
    return
end

[p, m] = size(sys);
n      = size(sys,'order');
l      = n / p;
l1     = l + 1;

%% Compute the observability matrix
O = zeros(l1 * p, n);
O(1:p, :) = sys.c;
for i = 2:l1
    O((i - 1) * p + 1:i * p, :) = O((i - 2) * p + 1:(i - 1) * p, :) * sys.a;
end

if (m > 0)
    P = null(O')';
    % Form the lower bock-triangular Toeplitz matrix T whose first 
    % block column is F = [ D; CB; CAB; ...; CA^{l-1}B ]
    F = [ sys.d; O(1:end-p, :) * sys.b ];
    T = zeros(l1 * p, l1 * m);
    for i = 1:l1
        T((i - 1) * p + 1:end, (i - 1) * m + 1: i * m) = F(1:(l1 + 1 - i) * p, :);
    end
    Q = P * T;
    R3 = [reshape(Q, p, m, l1) -reshape(P, p, p, l1)];
    % Reorder first exact then noisy
    ne = length(exct); nn = length(noisy);
    R  = [reshape(R3(:, exct,:), p, l1 * ne) reshape(R3(:, noisy, :), p, l1 * nn)]';
else % Output only (generically no exact variables)
    R = null(O');
end

%% Normalize with last block -I, in order to get X
x = -R(1:end - p, :) / R(end - p + 1:end, :);


%% INISTATE - compute initial condition for a trajectory
function xini = inistate(w, sys, T)

%% Define constants
[p, m] = size(sys.d);
n      = size(sys.a); n = n(1);
N      = size(w, 3);

if nargin < 3
    T = max(ceil(n / p), 2);
end

%% Define inputs and outputs
u  = w(1:T, 1:m, :);
y  = w(1:T, m + 1:end, :);
clear w

%% Compute the extended observability matrix
O = zeros(T * p, n);
O(1:p, :) = sys.c;
for t = 2:T
    O((t - 1) * p + 1:t * p, :) = O((t - 2) * p + 1:(t - 1) * p, :) * sys.a;
end

%% Compute initial consitions
xini = zeros(n, N);
sys.Ts = -1;
mo = (m > 0);
for k = 1:N
    if mo
        y0 = (y(:, :, k) - lsim(sys, u(:, :, k), 1:T))';
    else
        y0 = y(:, :, k)';
    end
    xini(:, k) = O \ y0(:);
end