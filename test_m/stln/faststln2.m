function [x,E,deltab,initflops,itflops,its,postflops]=faststln2(A,b,tol)

%This function is a fast STLN solver for
%toeplitz systems of the form Ax=b where
%[A b] is toeplitz. tol is a user defined
%tolerance indicating the degree to which
%rank deficiency is required or in other 
%words it precises the size of the smallest
%singular value of [A b]

noit=100;  %maximum number of iterations
[m,n]=size(A);
mn1=m+n-1;
T=[A b];

t=[T(1,n+1:-1:1).' ; T(2:m,1)];  

%Construction of e1
e1=[1 zeros(1,m-1)].';

%construction of F
F=zeros(m,mn1);
F(2:m,1:m-1)=eye([m-1,m-1]);

%initial estimate lambda
lambda=ones(m,1);
%lambda=zeros(m,1);
%initial estimate alpha
alpha=zeros(m+n-1,1);
beta=0;

%initial estimate x
initflops=0;
%flops(0)
x=A\b;
%init%flops=%flops;

if nargin < 4
tol = 1e-5;
end

%start iteration
%it%flops=0;
%post%flops=0;
for iii=1:noit
   %flops(0)
  %iii
  [rhs,rhat]=constrhs(alpha,beta,x,lambda,t,m,n);


 %keyboard
  %norm(rhat,2)

  if norm(rhat,2)<tol
    %Construction of E
    E=toeplitz(alpha(n:mn1),alpha(n:-1:1));
    %construction of deltab
    deltab=[beta;alpha(1:m-1)];
    %fprintf('Converged in %d iterations \n',iii)
    %post%flops=post%flops+%flops;
    break
  end

  
  alphapiu=alpha+t(2:m+n);

  %create generators
  g=genera(x,alphapiu,m,n);


  %[r,d]=fldstzb(g,m,n);
  [r,d]=fl2zbstab(g,m,n);


   %keyboard 
  v1=[m+n, 1:m+n-1,m+2*n+1:2*m+2*n,m+n+1:m+2*n];
  
  v2=[2:m+n, 1,2*m+n+1:2*m+2*n,m+n+1:2*m+n];
  b1=rhs(v1); 

  y=SolveL1b(x,m,n,b1);


  y(m+n+1:2*m+2*n)=SolL2b(r,m,n,y(m+n+1:2*m+2*n));


  y(m+n+1:2*m+n)=-1*y(m+n+1:2*m+n);
  y(m+n+1:2*m+2*n)=SolL2Tb(r,m,n,y(m+n+1:2*m+2*n));


  y=SolL1Tb(x,m,n,y);


  y=y(v2);y=y';
  %keyboard
  %update the variables
  alpha=alpha+y(1:m+n-1,1);
  beta=beta+y(m+n,1);
  x=x+y(m+n+1:m+2*n,1);
  lambda=lambda+y(m+2*n+1:2*m+2*n,1);
  %construction of  jacobian
  %E=toeplitz(alpha(n:m+n-1),alpha(n:-1:1));
  %J1=toeplitz([-x(n);zeros(m-1,1)],[-x(n:-1:1,1); zeros(m-1,1)]);
  %AA=[J1' ; [1 zeros(1,m-1)] ; -(A+E)'];  
  %gg=[alpha; beta; zeros(n,1)];  
  %lambda=-AA\gg;
%%flops,'yipi'   
  %check to see if converged
  %keyboard
  
  %it%flops=it%flops+%flops;
end
if iii==noit
    %Construction of E
    E=toeplitz(alpha(n:mn1),alpha(n:-1:1));
    %construction of deltab
    deltab=[beta;alpha(1:m-1)];
    %fprintf('Not Converged\n')
    %%%%post%flops=post%flops+%flops;

  %it%flops=it%flops/iii;
  %its=iii;
else
  %it%flops=it%flops/(iii-1);
  its=iii-1;
end

