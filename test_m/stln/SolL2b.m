
function [x]=SolL2b(r,m,n,b)
%
%solution of the triangular system r'x=b,
% where r is the Cholesky factor of M^T M
% Input : the matrix r, m,n
% output x:the solution of r'x=b;

mf=m+n;

x(1)=b(1)/r(1,1);
for i=2:n+1,
  sum=0;
  for j=1:i-1,
    sum=sum+r(j,i)*x(j);
  end
  x(i)=(b(i)-sum)/r(i,i);
end

for i=n+2:m,
  sum=0;
  for j=i-n:i-1,
    sum=sum+r(j,i)*x(j);
  end
  x(i)=(b(i)-sum)/r(i,i);
end


for i=m+1:n+m,
  sum=0;
  for j=1:i-1,
    sum=sum+r(j,i)*x(j);
  end
  x(i)=(b(i)-sum)/r(i,i);
end

x=x';