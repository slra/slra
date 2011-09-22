%This function calculates the structured Total Least Squares
%approximation of an extended Hankel matrix [A b].  
%Inputs:
%ynw: ynw=1 will minimize the Frobenius norm of (the approximating
%matrix-the original one). ynw=0 minimizes the difference between the
%signal and its approximation (signal can be read from first column
%and first row).
%[A b]: extended Hankel data matrix
%maxiter: maximum number of iterations
%tol: if convergence takes place, then norm(A*x-b)<tol
%startval:choose initialization method between 'ls' and 'htlsu'
%po: if po==1, detailed output during iteration is printed
%Output:
%[Atln Btln]: extended rank deficient Hankel matrix
%Xtln: parametervector
%fiter: number of iterations performed
%converged: if converged==1, STLN solution is found
%
%Author: Philippe Lemmerling
%Last modified: 16/08/1999

function [Xtln,iter,Ctln,fiter,converged] = stlnlserq(A,b,tol,maxiter,po)
	     
con=0;	    
[m,n]=size(A);
[m,d]=size(b);
mnd1=m+n+d-1;
mn1=m+n-1;  
converged=0;

ynw = 0;

%construction of D=[D_beta    0    ]
%                  [  0     D_alpha]
D=eye(mnd1,mnd1);
if ynw==1
  max_diaglength=min(n+d,m);
  D=D*max_diaglength;
  for i=1:max_diaglength
    D(i,i)=i;
    D(mnd1-i+1,mnd1-i+1)=i;
  end
end
D_beta=sqrt(D(1:d,1:d));
D_alpha=sqrt(D(d+1:mnd1,d+1:mnd1));

%construction of F
F=[];
for i=1:d
  Ftemp=zeros(m,d);
  Ftemp(max(1,m-i+1):m,d-i+1:min(d,d-i+m))=yey(min(m,i));
  F=[F;Ftemp];
end

%if strcmp(startval,'ls')
  x=tls(A,b);
  alpha=zeros(mn1,1);
  beta=0;
%elseif strcmp(startval,'htlsu')
  %put data into vector
%  hvect=[A(:,1);A(m,2:n).';b(m)];
%  [freq,damp,ampl,phas]=htlsu(hvect,hvect,n,1,1:m+n,ceil((m+n)/2)); 
%  htlsuvect=real(recoons(1:m+n,freq,damp,ampl,phas)).';
%  deltavect= htlsuvect-hvect;
%  alpha=deltavect(m+n-1:-1:1);
%  beta=deltavect(m+n);
%  Hhtlsu=hankel(htlsuvect(1:m),htlsuvect(m:m+n));
%  x=Hhtlsu(1:n,1:n)\Hhtlsu(1:n,n+1);
%end

fiter=x;

var=[beta;alpha;x];

%construction of GG
GG=[D_beta zeros(d,mn1) zeros(d,n*d);
   zeros(mn1,d) D_alpha zeros(mn1,n*d);
   zeros(n*d,mnd1+n*d)];

%start loop
iter=0;npr=1;
tol2=1e-2;
delta_var=100;
npr=tol+1;
while (iter<maxiter) %|(norm(delta_var)>tol2)
  iter=iter+1;
  if po==1
    fprintf('%d\n',iter)
    fprintf('Goal value %f\n',alpha.'*D_alpha*alpha+beta*D_beta*beta)
  end
  %construction of E matrix
  row=alpha(n:-1:1,1);
  col=alpha(mn1:-1:n,1);
  Etemp=hankel(col,row);
  E=zeros(m*d,n*d);
  for i=1:d
    E(1+(i-1)*m:i*m,1+(i-1)*n:i*n)=Etemp;
  end    
  	 
  %construction of matrix X
  X=[];  
  co=[]; 
  for i=1:d
    col=[zeros(m-1,1) ; x(n+(i-1)*n)];
    row=[ x(n*i:-1:1+(i-1)*n).' zeros(1,m-1) ];%keyboard
    Xtemp=hankel(col,row);
    corr=zeros(m,mn1);
    if (m-i)>0
      %corr=zeros(m,mn1);
      for jj=1:m-i
  	corr(jj,m-i-jj+1)=1;
      end
      Xtemp=Xtemp-corr;
    end  
    co=[co;corr];
    X=[X;Xtemp];
  end    
  %keyboard
  %construction of AA
  AA=[F' ; -X' ; -(A+E)'];
  	 
  %construction of vector gg
  gg=[D_beta*beta ;D_alpha*alpha;zeros(n*d,1)];
  	 
  %construction of vector res
  res=b+F*beta-A*x-X*alpha;
  %keyboard
  npr=norm(res);
  if po==1
    fprintf('Norm residue %g\n',npr)
  end
  
  if (npr<tol)&(delta_var<tol2)    
    converged=1;

    break;
  end    

  %S=AA';
  %S=S(m:-1:1,:);
  %s=res(m:-1:1);
  %[x]=lserq(GG,-GG*var,S,-s);
  %the above part was needed for the previous version of housen.m. That
  %version couldn't deal with x(1)==0 although this shouldn't be a problem
  %for applying the Householder transformation. I therefore had to permute 
  %the rows of the system of equations of the equality. The new version of
  %housen.m has no problem with x(1) being zero.
  %keyboard
  [x]=lserq(GG,-GG*var,AA',-res);
  
  %keyboard
  delta_var=x;	 
  var=var+delta_var;
  beta=var(1:d,1);
  alpha=var(d+1:d+mn1,1);
  x=var(d+mn1+1:d+mn1+n*d);
  fiter=[fiter x];
end

%keyboard
if converged ~= 1
%  fprintf('Terminated succesfully in %d iterations\n',iter);
%else
  fprintf('\nSTLN solution not found.\n');  
end

Atln=A+hankel(alpha(mn1:-1:n),alpha(n:-1:1));
Btln=reshape(co*alpha+F*beta+reshape(b,m*d,1),m,d);
Xtln=reshape(x,n,d);
%fiter=iter;

Ctln = [Atln Btln];