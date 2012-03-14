function [sysh, info, w, xini] = ident(w, m, l, opt)
% IDENT - Identification by structured low rank approximation
%
% [sysh, info, wh, xini] = ident(w, m, l, opt)
%
% Inputs:
% W   - set of trajectories, stored in an arry with dimensions
%                 #samples x #variables x #time series        
%       in case of equal number of samples or a cell array with  
%       #time series entries of dimension #samples x #variables
%       in case of different number of samples
% M   - input dimension
% L   - system lag
% OPT - options for the optimization algorithm:
%   OPT.EXCT - vector of indices for exact variables (default [])
%   OPT.SYS0 - initial approximation: an SS system with M inputs, 
%              P := size(W, 2) - M outputs, and order N : = L * P
%   OPT.DISP - level of display [off | iter | notify | final] (default notify)
%   OPT.MAXITER - maximum number of iterations (default 100)
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
    [T, nw, N] = size(w); % # samples x # variables x # time series
else
    [T{1}, nw] = size(w{1}); N = length(w);
    for i = 2:N
        T{i} = size(w{i}, 1);
        if size(w{i}, 2) ~= nw 
            error('All trajectories must have the same number of variables.')
        end
        if size(w{i}, 3) > 1
            error('The cell array entries of W must be 2D arrays.')
        end
    end
end
p  = nw - m;     % # of outputs
n  = p * l;      % order of the system
siso = (m == 1) & (p == 1);
l1   = l + 1; 

%% Check for parameters out of range
if nargin < 3
    error('Not enough input arguments. (W, M, and L should be given.)')
end
if ~iscell(w)
    if (~isreal(w)), error('Complex data not supported.'), end
else
    for i = 1:N
        if (~isreal(w{i})), error('Complex data not supported.'), end
    end
end
if (~isreal(m) | m < 0 | m > nw | (ceil(m) - m ~= 0))
    error('Invalid number of inputs M. (0 < M < size(W, 2) and integer.)')
end
if (~isreal(l) | l < 0 | (ceil(l) - l ~= 0))
    error('Invalid system lag L. (0 < L and integer.)')
end

%% Default options values
disp    = 'notify';
method  = 'l';
epsrel  = 1e-4;
epsabs  = 1e-4;
epsgrad = 1e-4;
maxiter = 1000;
sys0    = []; 
exct    = []; 

%% User specified options
OPTNAMES = {'disp', 'epsrel', 'epsabs', 'epsgrad', ...
            'maxiter', 'sys0', 'exct', 'method'};
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
opt = []; opt.maxiter = maxiter; opt.sys0 = sys0; opt.disp = disp; 
opt.epsrel = epsrel; opt.epsabs = epsabs; opt.epsgrad = epsgrad; 
opt.method = method;

%% Check for invalid option values
exct = unique(exct);
if any(1 > exct | nw < exct)
    error('An index for exact variable is out of range.')
end
if ~iscell(w)
    if (length(exct) * T > m * T + n)
        error('Too many exact variables. Generically there is no solution.')
    end
else
    for i = 1:N
        if (length(exct) * T{i} > m * T{i} + n)
            error('Too many exact variables. Generically there is no solution.')
        end
    end
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
if (N == 1) && (T <= n) 
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
if ~isa(opt.sys0, 'ss')
    warning('OPT.SYS0 not an SS object. Ignored.')
    sys.sys0 = [];
end
if ~isempty(opt.sys0)
    [pp, mm] = size(opt.sys0); 
    nn = size(opt.sys0, 'order');
    if (mm ~= m) | (pp ~= p) | (nn ~= n)
        warning('OPT.SYS0 invalid. Ignored.')
        sys.sys0 = [];
    end
end
% Modify the next line to change the default initial approximation 
if isempty(sys0), sys0 = n4sid(iddata(w(:, m + 1:end), w(:, 1:m)), l * p); end

%% Structure specification
ne = length(exct);
nn = nw - ne;
noisy = 1:nw; noisy(exct) = []; % indices of noisy variables
if ne == 0
    struct = [l1 nw]; 
    if ~iscell(w)
        if N > 1, struct.a = struct; struct.k = N; end
        par = vec(shiftdim(w, 1));
    else
        struct = kron(ones(N, 1), struct); par = [];
        for i = 1:N
            par = [par; 
                   vec(shiftdim(w{i}, 1))];
        end
    end
else % Define a block of exact variables
    struct = [l1 ne 1; l1 nn 0]; 
    if ~iscell(w)
        if N > 1, struct.a = struct; struct.k = N; end
        par = [vec(shiftdim(w(:, exct, :), 1));
               vec(shiftdim(w(:, noisy, :), 1))];
    else    
        struct = kron(ones(N, 1), struct); par = [];
        for i = 1:N
            par = [par;
                   vec(shiftdim(w(:, exct, :), 1));
                   vec(shiftdim(w(:, noisy, :), 1))];
        end
    end
end

%% Call the SLRA solver
x0 = sys2x(sys0, exct, noisy);
[x, info, cov, parh] = mex_slra(par, struct, l1 * m + n, x0, opt);
info.M = sqrt(info.fmin); info = rmfield(info, 'fmin');

if nargout > 2
    if ~iscell(w)
        if N == 1
            w(:, noisy) = reshape(parh(end - nn * T + 1:end), nn, T)';
        else % N > 1
            w(:, noisy, :) = shiftdim(reshape(parh(end - nn * T * N + 1:end), nn, N, T), 2);
        end
    else
        ind2 = 0; 
        for i = 1:N
            ind1 = ind2 + ne * T{i} + 1; 
            ind2 = ind2 + nw * T{i};
            w{i}(:, noisy) = reshape(parh(ind1:ind2), nn, T{i})';
        end
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
            tmp = w(1:l1, :, i)'; 
            xini(:, i) = O \ tmp(:);
        end
    end
end


%% VEC - vectorize a matrix
function a = vec(A), a = A(:);


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
if nargin < 3
    T = max(ceil(n / p), 2);
end

%% Define inputs and outputs
if ~iscell(w)
    N = size(w, 3);
    u = w(1:T, 1:m, :);
    y = w(1:T, m + 1:end, :);
else
    N = length(w);
    for i = 1:N
        u{i} = w{i}(1:T, 1:m, :);
        y{i} = w{i}(1:T, m + 1:end, :);
    end
end

%% Compute the extended observability matrix
O = zeros(T * p, n);
O(1:p, :) = sys.c;
for t = 2:T
    O((t - 1) * p + 1:t * p, :) = O((t - 2) * p + 1:(t - 1) * p, :) * sys.a;
end

%% Compute initial consitions
sys.Ts = -1; mo = (m > 0);
if ~iscell(w)
    xini = zeros(n, N); 
    for k = 1:N
        if mo
            y0 = (y(:, :, k) - lsim(sys, u(:, :, k), 1:T))';
        else
            y0 = y(:, :, k)';
        end
        xini(:, k) = O \ y0(:);
    end
else
    for k = 1:N
        xini{k} = zeros(n, N);
        if mo
            y0 = (y{k} - lsim(sys, u{i}, 1:T))';
        else
            y0 = y';
        end
        xini{k} = O \ y0(:);
    end
end
