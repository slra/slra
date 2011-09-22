%function [QT]=accq(R)
%
%This is a function specifically for the RQ factorization. In order to
%use this function, the previous call to rq has to be made with qunfact=0,
%otherwise R does not contain the Householder info. This function constructs
%in the most efficient way (i.e. starting with the "smallest" Householder
%transformations) the matrix QT, which is the same orthogonal transformation
%that has to be applied to the matrix in the original RQ factorization in
%order to obtain the R 
%
%Author: Philippe Lemmerling
%Last modified: 04/09/1999
%Reference: Golub & Van Loan p. 199, the so-called backward accumulation
%strategy 
function [QT]=accq(R)

[p,n]=size(R);
QT=eye([n,n]);
for i=1:p
   v=[1 ; R(i,1:n-p+i-1).'];
      %keyboard
   QT(1:n-p+i,1:n-p+i)=rowhouse(QT(1:n-p+i,1:n-p+i),v(n-p+i:-1:1,1));
end
