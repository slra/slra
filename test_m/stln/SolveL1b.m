function y = SolveL1b(x,m,n,b)

  x = [-1;x(n:-1:1)];
  y([1:m+n,2*m+n+1:2*m+2*n]) = b([1:m+n,2*m+n+1:2*m+2*n]);
  k = -1;
  for i = m+n+1:2*m+n,
    sum = 0;
    k = k+1;
    for j = 1:n+1,
      sum = sum + x(j)*y(k+j);
    end
    y(i) = b(i) - sum;
  end
  y = y';
