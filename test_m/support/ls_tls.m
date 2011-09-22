% LS-TLS - solves the Least Squares Total Least Squares problem  
% [A1 A2] X = B, where A1 is noise-free and [A2 B] is perturbed 
% by i.i.d. noise. 

function [X, Aest, Best, DelA, DelB] = wyLSTLS(A1, A2, B, r )

% ------------------------------------------------------------------
% get sizes of matrices
% ------------------------------------------------------------------

if isempty(B)
   error('B is an empty matrix!!!');
else
   [m d] = size(B);
end

if isempty(A1)
   m1 = m; n1 = 0;
else
   [m1 n1] = size(A1);
end

if isempty(A2)
   m2 = m; n2 = 0;
else
   [m2 n2] = size(A2);
end

if nargin < 4
   r = n2;
end

n = n1 + n2;  % convenient to have n

% ------------------------------------------------------------------
% checks
% ------------------------------------------------------------------
if m ~= m1 | m1 ~= m2 | m ~= m2
   error('Number of rows of matrices unequal'); 
else clear m1 m2;
end

% ------------------------------------------------------------------
% QR factorization
% ------------------------------------------------------------------

if n1 > 0
   [Q R] = wyHouseholder([A1 A2 B], n1);
   % partition into sub-matrices
   R11 = R(1:n1,1:n1);
   R12 = R(1:n1,n1+1:end);
   R22 = R(n1+1:end,n1+1:end);
   clear Q R;
   if n1 == n
      X = R11\R12; 
      Aest = A1;
      Best = A1*X;
      DelA = zeros(size(A1));
      DelB = B - Best;
      return;
   end
else
   R22 = [A1 A2 B];
end

% ------------------------------------------------------------------
% SVD
% ------------------------------------------------------------------

[U, S, V] = svd(R22);

% ------------------------------------------------------------------
% Rank determination: User specified
% ------------------------------------------------------------------

% ------------------------------------------------------------------
% Solution Space
% ------------------------------------------------------------------

if n1 > 0
   V22 = V(:,r+1:n2+d);
   V12 = R11\(-R12*V22);
   V2 = [V12;V22];
else
   V2 = V(:,r+1:n2+d);
end


% ------------------------------------------------------------------
% TLS solution
% ------------------------------------------------------------------

[Q R] = wyHouseholder(V2,n2-r);

Y = R(1:n,1:n2-r);
Z = R(1:n,n2-r+1:n2-r+d);
G = R(n+1:end,n2-r+1:end);

X = (-Z)/G;

DelAB = [A1 A2 B]*[Z; G]*[Z' G'];
DelA = DelAB(1:m, 1:n);
DelB = DelAB(1:m, n+1:n+d);
clear DelAB;

Aest = [A1 A2] - DelA;
Best = B - DelB;


function [Q,R] = wyHouseholder(A, n1);
% Compute Householder transformantion of a (mxn) matrix A till column n1
% and return a (mxm) orthogonal matrix Q, and a matrix R such that the 
% the first (n1xn1) sub-matrix is upper triangular.

% get the number of rows(m) and columns(n) of A
[m n] = size(A);

% set default: all columns worked on
if nargin < 2 n1 = n; end

% check
if m < n1, error('working on more columns than there are rows'); return; end

% create (mxm) Identity matrix, it will be used often
Im = eye(m);

% initialize R and Q
R = A;
Q = eye(m);

% loop over n1 columns
for k = 1:n1
   x = R(k:end,k);
   [v, beta] = gallery('house',x);
   v = [zeros(k-1,1); v];
   H = Im - beta*v*v';
   clear v beta;
   Q = Q*H;
   R = H*R;
end
clear k;

