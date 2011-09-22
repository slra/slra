function [l,d] = fl2zbstab(x,m,n)

  mn = m+n;
  m1 = m+1;
  mnmeno1 = mn-1;
  x1 = x(1,:);
  x2 = x(2,:);
  x3 = x(3,:);
  x4 = x(4,:);
  l(1,:) = x1;
  d(1) = -1;
  x1(2:mn) = x1(1:mn-1);
  x1(m+1) = 0;
  t3 = x3(m+1);
  t4 = x4(m+1);  
  t2 = x1(m+1);
  temp = 1;
  mmenon = m-n;

  for i = 2:mmenon
    i1 = i + 1;
    in = i + n;
    ni1 = n + i1;
    t2 = x1(m1);t4 = x4(m1);  
    [a,c,s] = giv(x1(i),x4(i));

    x4([i1:in, m1:mn]) = c*x4([i1:in, m1:mn]) - s*x1([i1:in, m1:mn]);
    
    x1(m1) = s*(t4-t3)/c;
    
    t3 = c*t3-s*x1(m1);     
    if i ~= mmenon,
      temp = temp*c;
      x4(ni1) = temp*x4(ni1);
    end
    s  = x2(i)/x1(i);  
    c1 = sqrt(x1(i)^2-x2(i)^2);
    c  = c1/x1(i);
    
    x1([i1:in,m:mn]) = (x1([i1:in,m:mn])-s*x2([i1:in,m:mn]))/c;
    x2([i1:in,m:mn]) = c*x2([i1:in,m:mn])-s*x1([i1:in,m:mn]);
    x1(i) = c1;
    
    l(i,i:mn) = x1(i:mn);
    d(i) = -1;
    x1(i1:mn) = x1(i:mnmeno1);
    x1(m1) = 0;
  end

  for i = m-n+1:m,
    i1 = i+1;
    t4 = x4(m1);  
    t2 = x1(m1);
    [a,c,s] = giv(x1(i),x4(i));

    x4([i1:mn]) = -s*x1([i1:mn])+c*x4([i1:mn]);
    x1(m1) = x1(m1)+s*(t4-t3)/c;
    t3 = c*t3-s*x1(m1);
    
    s  = x2(i)/x1(i);  
    c1 = sqrt(x1(i)^2-x2(i)^2);
    c  = c1/x1(i);
    
    x1([i1:mn]) = (x1([i1:mn])-s*x2([i1:mn]))/c;
    x2([i1:mn]) = -s*x1([i1:mn])+c*x2([i1:mn]);
    x1(i) = c1;

    l(i,i:mn) = x1(i:mn);
    d(i) = -1;
    x1(i1:mn) = x1(i:mnmeno1);
    x1(m1) = 0;
  end
  x3     = x4;
  x3(m1) = t3;

  for i = m1:mn,
    i1 = i+1;
    [a,c,s] = giv(x1(i),x4(i));
    t = c*x1(i:mn)+s*x4(i:mn);
    x4(i1:mn) = -s*x1(i1:mn)+c*x4(i1:mn);
    x1(i:mn) = t;

    [a,c,s] = giv(x2(i),x3(i));
    t = c*x2(i:mn)+s*x3(i:mn);
    x3(i1:mn) = -s*x2(i1:mn)+c*x3(i1:mn);
    x2(i:mn) = t;
    
    s  = x1(i)/x2(i);  
    c1 = sqrt(x2(i)^2-x1(i)^2);
    c  = c1/x2(i);
    
    x2([i1:mn]) = (x2([i1:mn])-s*x1([i1:mn]))/c;
    x1([i1:mn]) = -s*x2([i1:mn])+c*x1([i1:mn]);
    x2(i) = c1;

    l(i,i:mn) = x2(i:mn);
    d(i) = 1;
    x2(i1:mn) = x2(i:mnmeno1);
    x2(m1) = 0;
  end
