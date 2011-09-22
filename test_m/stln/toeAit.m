% file toeAit.m is the same as toeA.m but with VARIABLE number of iterations
% used !!!!! Here, niter is an output parameter that indicates the number of
% iterations performed by the TLN algorithm. The max. number of iterations
% is set to 50
%
% file toeAit.m for solving the d-dimensional TLN problem AX = B with A 
% Toeplitz and dimensions mxn, B unstructured of dimensions mxd.
% This program is for any Toeplitz matrix A where all the diagonals of A are 
% subject to error.
%
% The problem to be solved is formulated as follows:
% Given complex or real matrices A mxn, B mxd with A Toeplitz, find the
% TLN solution of AX = B by preserving the Toeplitz structure of A, i.e.
% find the BEST rank n structure preserving matrix approximation [A+E B+dB] 
% such that (A+E)X = B+dB is exactly solvable and (A+E) remains Toeplitz.
%
% the notation follows the paper by Rosen, Park, Glick except for r and the
% additional notes by Sabine Van Huffel
%
function[xtln,niter,E,dB]=toeAit(A,B,x,tol)

format short e
[m,n]=size(A);
[m,d]=size(B);
mn1=m+n-1; 
maxiter=50;
%
% compute d-dimensional TLN solution via iteration
% the perturbation E will have Toeplitz structure and r should be close 
% to zero
% initialize
% disp('program - tln algorithm, Toeplitz structure of A preserved'),
md=m*d; nd=n*d; pmax=min(m,n);
% set up the diagonal weighting matrix D that accounts for the repetition of
% elements of alpha (denoted by al here) in the matrix E 
%D=pmax*eye(mn1,mn1);
%for i=1:pmax,
% D(i,i)=i;
% D(m+n-i,m+n-i)=i;
%end;
%D=sqrt(D);
D = eye(mn1,mn1);
%
% disp('initial values :'),
% compute initial TLN solution by solving d ordinary LS problems Ax_i=b_i 

if (nargin < 3) | isempty(x)
x=zeros(n,d); 
for j=1:d,
    x(:,j)=A\B(:,j);
end;
end
if (nargin < 4)
  tol = 1e-10;
end

% disp('initial xtln: '), disp(x),
%
% initialize E, dB and al accordingly 
E=zeros(m,n);
al=zeros(mn1,1);
dB=A*x-B;

% now iterate
ncol=nd+mn1;
MM=zeros(md+mn1,ncol);
MM(md+1:md+mn1,1:mn1) = D;
colx=zeros(m,1); rowx=zeros(1,mn1);
del=ones(ncol); i=0;
%
% start iterations and stop when norm(del) is smaller than (nd+mn1)*10^(-15)
% or when the max. number of iterations is exceeded.
%
delx = inf;
nx0  = norm(x,'fro');
while (norm(delx,'fro')/nx0 > tol & i < maxiter)
    i=i+1;
%   disp('start of iteration'),disp(i),
    % Step 2(a): compute delal, delx

    % construct X
    X=[];
    for j=1:d,
        colx(1)=x(n,j);
        rowx(1:n)=x(n:-1:1,j);
        X0=toeplitz(colx,rowx);
        X=[X;X0];
    end;

    % vectorize dB to r
    r=reshape(dB,md,1); 
%   the following values are computed for test purpose only 
    rnorm=norm(r);
  % tlnerror IS the same as the total TLN correction norm([E dB],'fro')
  % tlnerror=sqrt(norm(r)^2 + (norm(D*al))^2);

  % construct variable part of MM (is matrix M of the paper)
    MM(1:md,1:mn1)=X; 
    for j=1:d,
        MM((j-1)*m+1:j*m,m+j*n:mn1+j*n) = A+E;
    end;
    del=-MM\[r ; D*al ];
  % disp('cond of M and norm(del):'), disp([cond(MM) norm(del)]),
    delal=del(1:mn1);
  % disp('norm delta alpha : '), disp(norm(delal)),
    delx=del(mn1+1:mn1+nd);
  % disp('norm delta x:  '), disp(norm(delx)),
    
    % Step 2(b): modify al and x 
    al=al+delal;
    x =x+reshape(delx,n,d);

    % Step 2(c): Construct new E from al and new dB 
    % form E
    cola=al(n:mn1);
    rowa=al(n:-1:1);
    E=toeplitz(cola,rowa);

    % form dB (negative residual value) 
    dB=(A+E)*x -B;
%   compute singular values of computed approximation [A+E,B+dB]
%   ss=svd([A+E, B+dB]); 
%   disp(' last min(m-n,d)+1 sing. values of [A+E,B+dB]='),disp(ss(n:min(m,n+d))'), 
end; % while loop 
niter=i;
% disp('final al: '), disp(al),
xtln=x;
% disp('final xtln: '), disp(xtln)

