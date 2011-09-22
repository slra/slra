 function[y]=vecty(M,m,n);
%
% construction of the vector y=[0;  M(2:2m+n-1,2:m+2n-1)^T M(2:2m+n-1,m+n)]
% input:M,m,n
% output:y
y(1)=0;
x=M(:,m+n);
 for i=2:n+1,
  sum=0;
  for j=2:i,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum+x(m+i-1);
 end  


 for i=n+2:m-1,
  sum=0;
  for j=i-n+1:i,
%  for j=2:i,
%    i,j,M(j,i),pause,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum+x(m+i-1);
 end  

 for i=m:m+n-1,
  sum=0;
  for j=i-n+1:m,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum+x(m+i-1);
 end  

 for i=m+n:m+2*n-1,
  sum=0;
  for j=2:m,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum;
 end  
