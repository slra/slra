function[x ]=Solverxb(r,m,n,b)
%
%solution of the triangular system rx=b,
% where r is the Cholesky factor of M^T M
%

mf=m+2*n-1;

x(mf)=b(mf)/r(mf,mf);
for i=mf-1:-1:mf-2*n,
  sum=0;
  for j=i+1:mf,
    sum=sum+r(i,j)*x(j);
  end
  x(i)=(b(i)-sum)/r(i,i);
end

for i=mf-2*n-1:-1:1,
  sum=0;
  for j=i+1:i+n-1,
    sum=sum+r(i,j)*x(j);
  end
  
  for j=m+n:mf,
    sum=sum+r(i,j)*x(j);
  end
 
  
  
  x(i)=(b(i)-sum)/r(i,i);
end
x=x';
