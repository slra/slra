 

function[c,s]=hyp(x,y);

if y==0,
  s=0.0;
  c=1.0;
elseif abs(y) < abs(x), 
  t=y/x;
  c=1.0/sqrt(1.0-t^2);
  s=c*t;
else abs(x) < abs(y),
  t=x/y;
  s=1.0/sqrt(1.0-t^2);
  c=s*t;
end
