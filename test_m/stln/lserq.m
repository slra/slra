%This function solves the equality constrained LS problem
%            min ||A*x-c||
%             x
%            s.t. B*x=d
%with A mxn, c mx1, B pxn, d px1
%using the generalized RQ factorization
%
%Author: Philippe Lemmerling
%Last modified: 16/08/1999
%Reference: LAPACK DGGLSE routine, see www.netlib.org for more info
function [x]=lserq(A,c,B,d)
[m,n]=size(A);
[p,dummy]=size(B);
M=p;N=(m-M)/2;
[dummy1,dummy2]=size(c);
[dummy3,dummy4]=size(d);

if (dummy~=n)|(dummy1~=m)|(dummy2~=1)|(dummy3~=p)|(dummy4~=1)
   error('Something is wrong with the dimensions of A, B, c or d. Type help lserq for more info')
end

%first part of generalized RQ (GRQ):
%RQ factorization of B
[T,Q]=rq(B,0);

%select triangular part of T 
%and solve for system by back substitution
x2=backsub(triu(T(:,n-p+1:n)),d);

%second part of GRQ:
%QR factorization of A*Q'
AQT=matxhous(A,T);
[Q,R]=qr(AQT);

QTc=Q'*c;
c1=QTc(1:n-p,1);
R11=R(1:n-p,1:n-p);
R12=R(1:n-p,n-p+1:n);

x1=backsub(R11,c1-R12*x2);

xT=matxhou2([x1'  x2'],T);
x=xT';

