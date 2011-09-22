function [r,c,s] = giv(a,b)       

  t = abs(a) + abs(b);
  if t == 0
    c = 1.0;
    s = 0.0;
  else
    ni = t*sqrt((a/t)^2 + (b/t)^2);
    c = a/ni;
    s = b/ni;
    r = ni;
  end

