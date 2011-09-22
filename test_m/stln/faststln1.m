function [x,iii,E,deltab]=faststln1(A,b,x,tol)
%This function is a fast STLN solver for
%toeplitz systems of the form Ax=b where
%A is toeplitz and b unstructured. tol is a user defined
%tolerance indicating the degree to which
%rank deficiency is required or in other 
%words it precises the size of the smallest
%singular value of [A b]

[pA,qA]=size(A);
slow=0; %setting this to 1 uses no optimization of efficiency at all
        %setting this to 0, you can choose between the displacement
        %rank approach (set method='disprank') or the rosen/park
        %optimization (set method='rospark')
method='disprank';
noit=150; %set the number of iterations in the STLN algorithm
%construction of D
D=eye([pA+qA-1,pA+qA-1]);

%initialize x and alpha
if nargin < 3 | isempty(x)
x=A\b;
end
alpha=zeros(pA+qA-1,1);

if nargin < 4
tol = 1e-5;
end

%start loop
for iii=1:noit 
  %construction of X
  X=[];
  for i=1:pA
    X=[X; zeros(1,i-1) x(qA:-1:1,1).' zeros(1,pA-i)];
  end
  %construction of E
  r=alpha(qA:-1:1);
  c=alpha(qA:pA+qA-1);
  E=toeplitz(c,r);
  
  %construction of M
  M=[X   A+E;
     D   zeros(pA+qA-1,qA)];

  %construction of r
  r=b-(A+E)*x;
 
  %construction of rhs
  rhs=[r;-D*alpha];

  if slow==1
    %flops(0)
    delta=M\rhs;
    %flops
  else
%    flops(0)
    y=[];
    y1=M(1,:);y1=y1';
    y2=[0,M(pA,1:pA+2*qA-2)];
    y2=y2';
    y(:,1)=y1;
    y(:,2)=y2;
    y1=M(pA+1,:);
    y(:,3)=y1';
    %MM=M'*M;
%    y1=M(2:2*pA+qA-1,2:pA+2*qA-1)'*M(2:2*pA+qA-1,pA+qA);
%    y1=[0;y1];
    y1=vecty(M,pA,qA);
    s=sqrt(y1(pA+qA));
    %s=sqrt(y1(pA+qA)+1);
    s1=s-1;
    y1=y1/s1;
    y1(pA+qA)=s;
    %keyboard
    y(:,4)=y1.';
    y1(pA+qA)=1;
    y(:,5)=y1.'; 
    %keyboard
%    tf0=flops;
%    fprintf('Generators: %d flops \n',tf0)      
    if strcmp(method,'disprank')
      %flops(0)
      R=terToep2(y,pA,qA);
      %[R]=FastP2(y,pA,qA);
      %tf1=flops;
      %fprintf('fast QR: %d flops \n',tf1)      
    elseif strcmp(method,'rospark')
%      flops(0)
      R=rrppgg(M,pA,qA);
%      tf1=flops;
%      fprintf('Rospark QR: 		%d flops \n',tf1)         
    elseif strcmp(method,'sqr')
%      flops(0)
      R=Rfact(M);
      R=R(1:pA+2*qA-1,:);
%      tf1=flops;
%      fprintf('Simple QR: %d flops \n',tf1)  
    elseif strcmp(method,'asqr')
%      flops(0)
      R=Rfact0(M);
      R=R(1:pA+2*qA-1,:);
%      tf1=flops;
%      fprintf('Advanced Simple QR: %d flops \n',tf1)       
    else
      error('You have initialized the string "method" with an invalid value')
    end
    %keyboard
    
    %the next line is an unoptimized lower triang. solver of me
    %temp2=lts(R.',M.'*rhs); 		%solve lower triang. system
    %flops(0)
    temp1=ProdScalMTx(M,pA,qA,rhs);
    %tf2=flops;
    %fprintf('rhs construct: %d flops \n',tf2)
    %flops(0)    
    %the next line is an optimized lower triang. solver by nicola
    [temp2]=Solvertxb(R,pA,qA,temp1);
    %tf3=flops;    
    %fprintf('LT system: %d flops \n',tf3)    
    %the next line is an unoptimized upper triang. solver of me
    %delta=uts(R,temp2); 		%solve upper triang. system
    %flops(0)    
    %the next line is an optimized upper triang. solver by nicola
    delta=Solverxb(R,pA,qA,temp2); 	%solve upper triang. system    
    %keyboard
    %tf4=flops;
    %fprintf('UT system: %d flops \n',tf4)    
  end
  dalpha=delta(1:pA+qA-1,1);
%  norm(delta)
  dx=delta(pA+qA:pA+2*qA-1,1);
  alpha=alpha+dalpha;
  x=x+dx;
  if norm(delta)<tol
    %fprintf('Convergence after %d iterations\n',iii)
    %keyboard
    break
  end
end

deltab=-r;