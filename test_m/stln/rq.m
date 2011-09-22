%function [R,Q]=rq(B,qunfact)
%
%This function calculates the RQ factorization of the pxn matrix B
%where p<=n:
%            B=R*Q=[0 Rs]*Q, 
%where Rs is upper triangular and square.
%If qunfact=1, then the orthogonal Q matrix is explicitly formed. If 
%qunfact~=1 then the successive householder transformations are stored
%below the diagonal of Rs implying that R is full. If Rs is needed it 
%can be obtained by selecting triu(R(:,n-p+1:n)).
%This version contains NO pivoting so if B is nearly rank deficient
%this program can not be used.
%
%Author: Philippe Lemmerling
%Last modified: 16/08/1999

function [R,Q]=rq(B,qunfact)
[p,n]=size(B);
if p>n
   error('This RQ factorization only deals with rectangular matrices having more columns than rows') 
end
R=B;%intialize triangular factor with B

if qunfact==1
   Q=eye([n,n]);
else
   Q=[];
end

for i=p:-1:1
   x=R(i,n-p+i:-1:1).';
   %Determine householder vector
   [top,v]=housen(x);
   %'hi'
   %keyboard
   
   %Apply Householder to row i of B
   R(i,n-p+i)=top;
   if qunfact~=1
     R(i,1:n-p+i-1)=v(2:n-p+i,1).';%below the diagonal the successive
                             %householder transformations are 
			     %stored in factored form (i.e. the 
			     %Householder vectors)
   else		     
     R(i,1:n-p+i-1)=zeros(1,n-p+i-1);	
     %if qunfact==1, R will be truly upperdiagonal, no Householder
     %info is stored subdiagonally
   end
  
   %Apply Householder to rows i-1:-1:1
   if i>1
      %R(1:i-1,n-p+i:-1:1)=colhouse(R(1:i-1,n-p+i:-1:1),v);
      R(1:i-1,1:n-p+i)=colhouse(R(1:i-1,1:n-p+i),v(n-p+i:-1:1));
      %R(:,n-p+i:-1:1)=colhouse(R(:,n-p+i:-1:1),v);
      %keyboard
   end
   
   if qunfact==1
      %Q(:,n-p+i:-1:1)=colhouse(Q(:,n-p+i:-1:1),v);
      Q(:,1:n-p+i)=colhouse(Q(:,1:n-p+i),v(n-p+i:-1:1));
   end
end
if qunfact==1
  Q=Q';
end

      
   
