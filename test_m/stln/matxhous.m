%function [AQT]=matxhous(A,R)
%
%This is a function specifically for the RQ factorization. In order to
%use this function, the previous call to rq has to be made with qunfact=0,
%otherwise R does not contain the Householder info. This function applies
%the same sequence of transformations applied to B in rq.m, but now the 
%transformations are applied to A. Therefore A needs to have the same 
%number of columns as B (since the Householder transformation in rq.m is a 
%post-multiplication in matrix terms). 
%
%Author: Philippe Lemmerling
%Last modified: 16/08/1999

function [AQT]=matxhous(A,R)
[m,dummy]=size(A);
[p,n]=size(R);
if dummy~=n
   error('The Householder transformation you try to apply is not compatible with A. Wrong number of columns! Type help matxhous.')
end

AQT=A;
for i=p:-1:1 %retrieve the consecutive Householder vectors
   v=[ 1 ; R(i,1:n-p+i-1).'];
   AQT(:,1:n-p+i)=colhouse(AQT(:,1:n-p+i),v(n-p+i:-1:1,1));
end
