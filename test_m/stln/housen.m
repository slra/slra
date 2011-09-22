%function [top,v]=housen(x)
%
%This function computes the householder vector v s.t.
%(I-2*v*v'/(v'*v))*x=[top 0 0 ... 0 0].'
%                         \_________/
%                         length(v)-1
%AND
%v(1) is normalized to 1 
%
%Author: Philippe Lemmerling
%Last modified: 16/08/1999
%Reference: GVL, p.196

function [top,v]=housen(x)
n=length(x);
mu=norm(x,2);
v=x;
if mu~=0
   beta=x(1)+(sign(x(1))+(x(1)==0))*mu;
   v(2:n,1)=v(2:n,1)/beta;
else
   error('Householder reflection applied to 0 vector!')
end
v(1)=1;
top=-(sign(x(1))+(x(1)==0))*mu;

