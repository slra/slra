function [x]=SolverTxb(r,m,n,b)
%
%solution of the triangular system r'x=b,
% where r is the Cholesky factor of M^T M
% Input : the matrix r, m,n
% output x:the solution of r'x=b;

mf=m+2*n-1;

x(1)=b(1)/r(1,1);
for i=2:n,
  sum=0;
  for j=1:i-1,
    sum=sum+r(j,i)*x(j);
  end
  x(i)=(b(i)-sum)/r(i,i);
end

for i=n+1:n+m-1,
  sum=0;
  for j=i-n+1:i-1,
    sum=sum+r(j,i)*x(j);
  end
  x(i)=(b(i)-sum)/r(i,i);
end


for i=m+n:2*n+m-1,
  sum=0;
  for j=1:i-1,
    sum=sum+r(j,i)*x(j);
  end
  x(i)=(b(i)-sum)/r(i,i);
end

x=x';