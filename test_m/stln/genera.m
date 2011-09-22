% computes the generator matrix for M
% input: the vectors x, \alpha and the dimensions m, n
% output: g, the generator matrix wrt the block shift matrix.

function[g]=genera(x,alpha,m,n)

  vz1 = ones(m-1,1);
  vz2 = ones(n-1,1);
  VZ1 = diag(vz1,-1);
  VZ2 = diag(vz2,-1);
  VZ3 = zeros(m,n);
  v   = zeros(m+n,1);

  x(n+1) = -1;
  y = x/norm(x,2);
  l = -1;
  for i = n+1:-1:1,
    sum = 0;
    l = l+1;
    v(l+1) = -x(i+l:-1:1+l)'*y(i:-1:1); 
  end
  v(m+1:m+n) = alpha(n:-1:1);

  xt = (v(1));
  x1(1) = xt;
  x2(1) = 0;
  for i = 2:n+1,
    x1(i) = v(i);
    x2(i) = x1(i);
  end

  for i = m+1:m+n,
    x1(i) = -v(i)/xt;
    x2(i) = x1(i);
  end

  v = zeros(m+n,1);
  v(1:m) = alpha(n:m+n-1);
  v(1) = 0;

  xt = .5;
  x3(m+1) = .5;
  x4(m+1) = -.5;
  for i = 2:m+n,
    if i ~= m+1,
      x3(i) = v(i);
      x4(i) = x3(i);
    end
  end
  g = [x1;x2;x3;x4];
