% Fast STLN solver for Toeplitz systems

function [x,its,E,deltab] = stls_n(A,b,x,tol)

  noit  = 500;  %maximum number of iterations
  [m,n] = size(A);
  mn1   = m+n-1;
  T     = [A b];

  t  = [T(1,n+1:-1:1)' ; T(2:m,1)];  
  e1 = [1 zeros(1,m-1)]';
  F  = zeros(m,mn1);
  F(2:m,1:m-1) = eye([m-1,m-1]);

  %initial estimate lambda
  lambda = ones(m,1);

  %initial estimate alpha
  alpha = zeros(m+n-1,1);
  beta  = 0;

  %initial estimate x
  if (nargin < 3) | (isempty(x))
    x = tls(A,b);
  end
  if (nargin < 4)
    tol = 1e-3;
  end
  
  nx0 = norm(x);
  x_  = zeros(n,1);
  %start iteration
  for iii = 1:noit
    [rhs,rhat] = constrhs(alpha,beta,x,lambda,t,m,n);
    if norm(x-x_)/nx0 < tol
      %Construction of E
      E = toeplitz(alpha(n:mn1),alpha(n:-1:1));
      %construction of deltab
      deltab = [beta;alpha(1:m-1)];
      %fprintf('Converged in %d iterations \n',iii)
      break
    end
    alphapiu = alpha + t(2:m+n);
    %create generators
    g = genera(x,alphapiu,m,n);
    %[r,d] = fldstzb(g,m,n);
    [r,d] = fl2zbstab(g,m,n);
    v1 = [m+n, 1:m+n-1,m+2*n+1:2*m+2*n,m+n+1:m+2*n];
    v2 = [2:m+n, 1,2*m+n+1:2*m+2*n,m+n+1:2*m+n];
    b1 = rhs(v1); 
    y  = SolveL1b(x,m,n,b1);
    y(m+n+1:2*m+2*n) = SolL2b(r,m,n,y(m+n+1:2*m+2*n));
    y(m+n+1:2*m+n)   = -1*y(m+n+1:2*m+n);
    y(m+n+1:2*m+2*n) = SolL2Tb(r,m,n,y(m+n+1:2*m+2*n));
    y = SolL1Tb(x,m,n,y);
    y = y(v2); 
    y = y';
    %update the variables
    alpha  = alpha+y(1:m+n-1,1);
    beta   = beta+y(m+n,1);
    x_     = x;
    x      = x+y(m+n+1:m+2*n,1);
    lambda = lambda+y(m+2*n+1:2*m+2*n,1);
    %construction of  jacobian
  end

  if iii == noit
    %Construction of E
    E = toeplitz(alpha(n:mn1),alpha(n:-1:1));
    %construction of deltab
    deltab = [beta;alpha(1:m-1)];
    fprintf('Not Converged\n')
    its = inf;
  else
    its = iii - 1;
  end

