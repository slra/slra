           
function[r]=terToep2(y,m,n);
mf=max(size(y));
mf1=mf-1;
y1=y(:,1)';
y2=y(2:mf,2)';
y3=y(:,3)';
y4=y(2:mf,4)';
y5=y(2:mf,5)';
%
% calcolo della matrice triangolare superiore r
% computation of the upper triangular matrix r

mn=m+n;
mn1=mn-1;
% 
a4=y(mn,4);a5=y(mn,5);vx=1;
%
% Computation of the firs row of r
%
[y1(1),c,s]=giv(y1(1),1);
x(2:n)=c*y1(2:n);
x(mn:mf)=c*y1(mn:mf);
y3(1:n-1)=-s*y1(2:n);
y3(mn-1:mf-1)=-s*y1(mn:mf);
y1(2:n)=x(2:n);
y1(mn:mf)=x(mn:mf);
r(1,1:mf)=y1;


for i=1:m-1,
ni1=n+i-1;
i1=i+1;  
 % updating  (x_1, x_3) 
  [x(i),c,s]=giv(y1(i),y3(i));
  x(i:ni1)=c*y1(i:ni1)+s*y3(i:ni1);
  x(mn1:mf1)=c*y1(mn1:mf1)+s*y3(mn1:mf1);
  y3(i1:ni1)=-s*y1(i1:ni1)+c*y3(i1:ni1);
  y3(mn1:mf1)=-s*y1(mn1:mf1)+c*y3(mn1:mf1);
  y1(i:ni1)=x(i:ni1);
  y1(mn1:mf1)=x(mn1:mf1);
 
  [a,c,s]=giv(y1(i),y4(i));
  a=y1(mn1)+s/c*(a4-a5);
  au=y1(mn1); 
 
  y4(i1:ni1)=-s*y1(i1:ni1)+c*y4(i1:ni1);
  y4(mn1:mf1)=-s*y1(mn1:mf1)+c*y4(mn1:mf1);
  y1(mn1)=a;
  
  a=a4;
  a5=-s*au+(a5-s^2*a4)/c;
  a4=y4(mn1);
  
  vx=c*vx;
  if i ~= m-1,
    y4(n+i)=vx*y4(n+i);
  end 

 
  r(i1,i1:mf)=y1(i:mf1);
  y1(i1:mf)=y1(i:mf1);
%'hi'
%keyboard 
end

for i=m:mn-2,
    i1=i+1; 
   [x(i),c,s]=giv(y1(i),y3(i));
  x(i1:mf1)=c*y1(i1:mf1)+s*y3(i1:mf1);
 
  y3(i1:mf1)=-s*y1(i1:mf1)+c*y3(i1:mf1);
 
  y1(i:mf1)=x(i:mf1);
 

  [a,c,s]=giv(y1(i),y4(i));
  a=y1(mn1)+s/c*(a4-a5);
  au=y1(mn1); 
 
  y4(i1:mf1)=-s*y1(i1:mf1)+c*y4(i1:mf1);
  
  y1(mn1)=a;
 
  a=a4;
  a5=-s*au+(a5-s^2*a4)/c;
  a4=y4(mn1);

   [c,s]=hyp(y1(i),y2(i));
   x=c*y1(i:mf1)-s*y2(i:mf1);
  y2(i1:mf1)=-s*y1(i1:mf1)+c*y2(i1:mf1);
  y1(i1:mf)=x;
 
  r(i1,i1:mf)=y1(i1:mf);
 
end
y5(mn1)=a5;
y5(mn:mf1)=y4(mn:mf1);

for i=mn1:mf1,
 i1=i+1;
  [a,c,s]=giv(y1(i),y3(i));
  x=c*y1(i:mf1)+s*y3(i:mf1);
  y3(i1:mf1)=c*y3(i1:mf1)-s*y1(i1:mf1);
  y1(i:mf1)=x;

  [a,c,s]=giv(y1(i),y4(i));
  x=c*y1(i:mf1)+s*y4(i:mf1);
  y4(i1:mf1)=c*y4(i1:mf1)-s*y1(i1:mf1);
  y1(i:mf1)=x;
 
    [a,c,s]=giv(y5(i),y2(i));
     x=c*y5(i:mf1)+s*y2(i:mf1);
     y2(i1:mf1)=c*y2(i1:mf1)-s*y5(i1:mf1);
     y5(i:mf1)=x;
 
  [c,s]=hyp(y1(i),y5(i));
  x=c*y1(i:mf1)-s*y5(i:mf1);
  y5(i1:mf1)=c*y5(i1:mf1)-s*y1(i1:mf1);
  y1(i1:mf)=x;
  
  r(i1,i1:mf)=x;

end
