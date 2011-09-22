%deconvB.m *******************************************************************
diary tlntest.out
diary on
format short e
% set up initial data for a deconvolution set of equations A and B
% A is an mxn matrix of white noise input samples (of mean 0 and var.20)
% the solution H is some finite impulse response (see fig.7.1.b. in doctoral
%                          thesis but I multiplied each value with 10)
%
%noisestd=[ 1.0e-10 1.0e-8 1.0e-6 1.0e-4 1.0e-2 5.0e-2 1.0e-1 5.0e-1];
% here, we consider noise std. instead of the noise variances considered in 
% the test example doct. thesis fig.7.1.b
%noisestd=[ 1.0e-10 1.0e-8 1.0e-6 1.0e-4 5.0e-4 1.0e-3 5.0e-3 1.0e-2 5.0e-2 1.0e-1 2.0e-1 3.0e-1 ];
noisestd=[ 1.0e-10 1.0e-8 1.0e-6 1.0e-04 1.0e-02];
ncases=length(noisestd);
m=30;
n=20;
d=3;
niter=10;
ntest=1;
pmax=min(m,n+d);
disp(' deconvolution example with A Toeplitz and non-zero upper triangular part')
disp('-----------------------------------------------------------------------')
%inputex=sqrt(20.0)*randn(m+n-1,1);
inputex=randn(m+n-1,1);
Aex=toeplitz(inputex(n:n+m-1),inputex(n:-1:1));
xex=zeros(n,d);
xex(:,1)=[0.19;0.33;0.44;0.54;0.59;0.62;0.63;0.61;0.58;0.56;0.53;0.50;0.485;0.46;0.40;0.34;0.18;0.10;0.02;0.0];
xex(:,2)=[0.29;0.43;0.64;0.74;0.89;0.92;0.93;0.91;0.78;0.76;0.73;0.60;0.585;0.56;0.50;0.44;0.28;0.10;0.02;0.0];
xex(:,3)=[0.09;0.23;0.34;0.44;0.49;0.52;0.58;0.56;0.53;0.51;0.48;0.45;0.40;0.36;0.30;0.24;0.08;0.010;0.002;0.0];
xex=10*xex;
nex=norm(xex,'fro');
disp('norm(xex)='),nex
disp('xex='),disp(xex),
[u,ss,v]=svd(Aex,0);
disp(' sing.values of Aex= '),diag(ss)'
disp(' xex transpose v= '),xex'*v
% compute output vector Y
Bex=Aex*xex;
xerr0=zeros(ncases,1);
xerr1=zeros(ncases,1);
xerr2=zeros(ncases,1);
lscorr=zeros(ncases,1); 
tlscorr=zeros(ncases,1); 
tlncorr=zeros(ncases,1); 
truecorr=zeros(ncases,1); 

disp(' ')
disp('m,n,d,ntest,niter='), disp ([m,n,d,ntest,niter]),
disp(' ')
disp('average solution error obtained by each method') 
disp('   noisestd    LS sol.err.  TLS sol.err.  TLN sol.err.  TLN noise s.v.')
disp(' --------------------------------------------------------------------- ')

for l=1:ncases;
noisstd=noisestd(l);
tlnsv=0.0;
for I=1:ntest;
%add noise to the exact data input , Y
input=inputex + noisstd*randn(m+n-1,1);
B=Bex+noisstd*randn(m,d);
cola=input(n:n+m-1);
rowa=input(n:-1:1);
A=toeplitz(cola,rowa);
corr=norm([A-Aex B-Bex],'fro');
truecorr(l)=truecorr(l)+corr;
%
% compute LS solution
%
xls=A\B;
corr=norm(A*xls-B,'fro');
% disp(' total LS correction=norm(Axls-B):'), corr,
lscorr(l)=lscorr(l)+corr;
%
% compute TLS solution via SVD
%
% the perturbation will not have Toeplitz structure
% disp('program - TLS algorithm - solution via SVD'),
[u,s,v]=svd([A B],0);
% this is a simple program where we assume that the TLS solution exists
% xtls is the tls solution
xtls=-v(1:n,n+1:n+d)*inv(v(n+1:n+d,n+1:n+d));
% disp('xtls: '), disp(xtls),
% Etls and dBtls are the perturbation on A and B, respectively
% (A+Etls)*xtls=B+dBtls (mathematically)
Etls=-u(:,n+1:pmax)*s(n+1:pmax,n+1:pmax)*v(1:n,n+1:pmax)';
dBtls=-u(:,n+1:pmax)*s(n+1:pmax,n+1:pmax)*v(n+1:n+d,n+1:pmax)';
% compute the TLS residue (not the true residue but here equal to dBtls)
% rtls=(A+Etls)*xtls-B; % disp('rtls: '), disp(rtls),
corr=norm(s(n+1:pmax,n+1:pmax),'fro');   %corr always=norm([Etls dBtls],'fro')
% disp('total TLS correction=norm([Etls dBtls]): ' ), disp(corr),
tlscorr(l)=tlscorr(l)+corr;
%
% The following x should be the same as xtls
% it is computed only for the purpose of verifying the computation
% xtls2=(A+Etls)\(B+dBtls);
% disp(' check: xtls-xtls2= '), norm(xtls-xtls2,'fro'),
%
% compute TLN solution
%
[xtln,E,dB]=toeA(A,B,niter);
ss=svd([A+E B+dB]);
tlnsv=tlnsv+norm(ss(n+1:pmax));
% disp(' last min(m-n,d)+1 sing.values:'),ss(n:pmax)',
corr=norm([E dB],'fro');
% disp('total TLN correction = norm([Etln,dBtln]): ' ), disp(corr),
tlncorr(l)=tlncorr(l)+corr;
%
% compare accuracy ls, tls and tln solutions
%
xerrls=norm(xls-xex,'fro')/nex;
xerr0(l)=xerr0(l)+xerrls;
xerrtls=norm(xtls-xex,'fro')/nex;
xerr1(l)=xerr1(l)+xerrtls;
xerrtln=norm(xtln-xex,'fro')/nex;
xerr2(l)=xerr2(l)+xerrtln;
end %loop I=1:ntest
%
% compute average relative errors
%
xerr0(l)=xerr0(l)/ntest;
xerr1(l)=xerr1(l)/ntest;
xerr2(l)=xerr2(l)/ntest;
tlnsv=tlnsv/ntest;
%
% display relative errors of tls solution and tln solution for
% each noise std considered
 fprintf(' %e  %e  %e %e  %e\n',noisestd(l),xerr0(l),xerr1(l),xerr2(l),tlnsv);
end %loop l=1:ncases
%
% display table with the total corrections applied to [A B] by each method 
%
lscorr=lscorr/ntest;
tlscorr=tlscorr/ntest;
tlncorr=tlncorr/ntest;
truecorr=truecorr/ntest;
disp(' ')
disp(' average norm of total corrections applied to [A B]') 
disp('   noisestd     LS corr.     TLS corr.     TLN corr.     true corr.  ')
disp(' --------------------------------------------------------------------- ')
for l=1:ncases
 fprintf(' %e  %e  %e %e  %e\n',noisestd(l),lscorr(l),tlscorr(l),tlncorr(l),truecorr(l));
end % loop l=1:ncases;

diary off
