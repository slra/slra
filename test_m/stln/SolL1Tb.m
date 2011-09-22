function y = SolL1Tb(x,m,n,b)

  x(n+1) = -1;
  y(m+n+1:2*m+2*n) = b(m+n+1:2*m+2*n);
  k = 2*m + n;
  l = m + n + 1;
  for i  =  m+n:-1:m+1,
    k = k-1;
    sum = 0;
    for j = 1:l-i,
      sum = sum + x(j)*y(k+j);
    end
    y(i) = b(i) - sum;
  end
  for i = m:-1:n+1,
    k = k-1;
    sum = 0;
    for j = 1:n+1,
      sum = sum + x(j)*y(k+j);
    end
    y(i) = b(i) - sum;
  end
  l = 0;
  for i = n:-1:1,
    sum = 0;
    l = l + 1;
    for j = 1:i,
      sum = sum+x(l+j)*y(k+j);
    end
    y(i) = b(i) - sum;
  end


