% file toeAB.m for solving the d-dimensional TLN problem with [A B] Toeplitz 
% this program is for any Toeplitz matrix [A B] where all
% the diagonals of [A B] are subject to error
%
% The problem to be solved is formulated as follows:
% Given complex or real matrices A mxn, B mxd with [A B] Toeplitz, find the
% TLN solution of AX = B by preserving the Toeplitz structure of [A B], i.e.
% find the BEST rank n structure preserving Toeplitz matrix approximation
% [A+E B+dB] such that (A+E)X = B+dB is exactly solvable.
%
% niter is the number of iterations performed by the tln algorithm
% the notation follows the paper by Rosen, Park, Glick except for r and the
% additional notes by Sabine Van Huffel
%
function [xtln,iter,E,dB,D] = toeAB(A,B,x,tol)

max_iter = 500;

format short e
[m,n]=size(A);
[m,d]=size(B);
mnd1=m+n+d-1; 
% W is the weighting factor used for solving the equivalent unconstrained LS
% problem. A large W reduces \hat r to zero in a very fast way but also
% makes MM ill-conditioned. 
W=10^5;
%
% compute d-dimensional TLN solution via iteration
% the perturbation [E dB] will have Toeplitz structure and r should be close 
% to zero
% initialize
md=m*d; nd=n*d; pmax=min(m,n+d);
% set up the diagonal weighting matrix D that accounts for the repetition of
% elements of alpha and beta (denoted by al and be here) in the matrix E and dB
%D=pmax*eye(mnd1,mnd1);
%for i=1:pmax,
% D(i,i)=i;
% D(m+n+d-i,m+n+d-i)=i;
%end;
%D = sqrt(D);
D=eye(mnd1); %sqrt(D);
% set up the matrix \hat P0 and denoted here by P0
P0=zeros(md,m+n-1);
for j=1:min(m-1,d),
    P0((j-1)*m+j+1:j*m,1:m-j)=eye(m-j,m-j);
end;
% set up the matrix \hat e0 and denoted here by e0
e0=zeros(md,d);
for j=1:d,
    minmj=min(m,j);
    e0((j-1)*m+1:(j-1)*m+minmj,d-j+1:d-j+minmj)=eye(minmj,minmj);
end;
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
  tol = 1e-3;
end

% initialize E, dB, r and al, be accordingly 
E=zeros(m,n);
dB=A*x-B;
r = reshape(dB,md,1);
al=zeros(m+n-1,1);
be=zeros(d,1); 

% now iterate
% an integer value for iter should be given 
MM=zeros(md+mnd1,nd+mnd1);
MM(1:md,1:d) = -W*e0; MM(md+1:md+mnd1,1:mnd1) = D;
colx=zeros(m,1); rowx=zeros(1,m+n-1);
x_  = zeros(n,d);
nx0 = norm(x,'fro');
%
% start iterations and stop when the given number of iterations is exceeded.
%
iter = 0;
delx = inf;
while (norm(delx,'fro')/nx0 > tol) & (iter < max_iter)
  iter = iter + 1;
%      disp('start of iteration'),disp(i),
    
    % Step 2(a): compute delal, delbe, delx

    % construct X
    X=[];
    for j=1:d,
        colx(1)=x(n,j);
        rowx(1:n)=x(n:-1:1,j);
        X0=toeplitz(colx,rowx);
        X=[X;X0];
    end;

  % the following values are computed for test purpose only 
  %  rnorm=norm(r);
  % tlnerror is NOT the same as the total TLN correction norm([E dB],'fro')!!
  % tlnerror=sqrt(norm(r)^2 + (norm(D(d+1:mnd1,d+1:mnd1)*al))^2);
  % disp('norm of r and tlnerror: '), disp([rnorm tlnerror]),

    % construct variable part of MM (is matrix M of the paper)
    MM(1:md,d+1:mnd1) = W*(X-P0); 
    for j=1:d,
        MM((j-1)*m+1:j*m,m+d+j*n:mnd1+j*n) = W*(A+E);
    end;
    del=-MM\[W*r ; D*[be;al] ];
%   disp('cond of M and norm(del):'), disp([cond(MM) norm(del)]),
    delal=del(d+1:mnd1);
%   disp('norm delta alpha : '), disp(norm(delal)),
    delbe=del(1:d);
%   disp('norm delta beta: '), disp(norm(delbe)),
    delx=del(mnd1+1:mnd1+nd);
%   disp('norm delta x:  '), disp(norm(delx)),
    
    % Step 2(b): modify al, be and x 
    al = al+delal;
    be = be+delbe;
    x_ = x;
    x  = x+reshape(delx,n,d);

    % Step 2(c): Construct new E from al and new dB from be and compute \hat r
    % form E
    cola = al(n:n+m-1);
    rowa = al(n:-1:1);
    E = toeplitz(cola,rowa);

    % form dB 
    dB = toeplitz([be(d);al(1:m-1)],be(d:-1:1));
%   compute singular values of computed approximation [A+E,B+dB]
%   ss=svd([A+E, B+dB]); 
%   disp(' last d+1 sing. values of [A+E,B+dB]='), disp(ss(n:n+d)'), 

    % construct (A+E)*x - (B+dB) and vectorize to r
    r = reshape((A+E)*x-B-dB,md,1);
end; % loop i=1:niter
% disp('final al: '), disp(al),
xtln=x;
% disp('final xtln: '), disp(xtln)

