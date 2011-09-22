  function[y]=ProdScalMTx(M,m,n,x);
  % Computatin of
  % y=M^T x
  % where M is the matrix [X,E;I,0]
  %
  
  
  for i=1:n,
  sum=0;
  for j=1:i,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum+x(m+i);
 end  


 for i=n+1:m-1,
  sum=0;
  for j=i-n+1:i,
%  for j=2:i,
%    i,j,M(j,i),pause,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum+x(m+i);
 end  

 for i=m:m+n-1,
  sum=0;
  for j=i-n+1:m,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum+x(m+i);
 end  

 for i=m+n:m+2*n-1,
  sum=0;
  for j=1:m,
    sum=sum+M(j,i)*x(j);
  end
  y(i)=sum;
 end  
y=y';