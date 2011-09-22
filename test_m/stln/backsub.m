%function [x]=backsub(R,b)
%
%backsubstitution for solving the linear system 
%            Rx=b
%in which R is square upper triangular
%
%Author: Philippe Lemmerling
%Last modified: 16/08/1999
function [x]=backsub(R,b)
[n,dummy]=size(R);
for i=n:-1:1
   temp=b(i);
   for j=i+1:n
      temp=temp-R(i,j)*x(j);
   end
   x(i,1)=temp/R(i,i);
end
