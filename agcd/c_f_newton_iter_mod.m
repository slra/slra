% This is a modified file from the software 
% Fastgcd by Paola Boito, paola.boito AT unilim.fr
%
% The original software should be downloaded from
% http://www.unilim.fr/pages_perso/paola.boito/software.html
%
% This file contains an extract of the file fast_gcd_wls.m,
% starting from the function c_f_newton_iter.
% The modified function was renamed as c_f_newton_iter_mod,
% to adapt it for speed comparison and per-iteration comparison.
%
% The modifications were made by Konstantin Usevich konstantin.usevich AT statmod.ru
%
% -------------------------------------------------
% function c_f_newton_iter_mod
% -------------------------------------------------
%
% Computes fast iterative refinement using Newton's method.
%
% [z,res]=c_f_newton_iter(f,g,z,q1,q2)
%
% f, g -> polynomials whose GCD is sought,
% z -> tentative GCD (to be refined),
% q1, q2 -> tentative cofactors (to be refined).
%
function [z,q1,q2,res, num_iter]=c_f_newton_iter_mod(f,g,z,q1,q2,optflag, max_num_iter, tol);

K=norm(z)^2+norm(q1)^2+norm(q2)^2;
% set tolerance for nearness and maximum number of iterations
%tol=1e-15;
%max_num_iter=20;
grado=length(z)-1;
N=length(f)-1;
M=length(g)-1;
f_nearness = [norm(conv(q1,z)-f), norm(conv(q2,z)-g)];
allow=1;
num_iter=0;
while allow==1
    num_iter=num_iter+1;
    altro_Z=c_iterfast(f,g,q1,q2,z,K,optflag);
    altro_z=fliplr((altro_Z(1:grado+1)).'); % new gcd
    altro_q1=fliplr((altro_Z(grado+2:2+N)).');
    altro_q2=fliplr((altro_Z(N+3:N+3+M-grado)).'); % new cofactors
    s_nearness=[norm(conv(altro_q1,altro_z)-f), norm(conv(altro_q2,altro_z)-g)];
    if (abs(norm(s_nearness)^.2 - norm(f_nearness)^.2)<tol)||(num_iter>max_num_iter)
%    if (norm(s_nearness)<tol)||(num_iter>max_num_iter)
        allow=0;
        z=altro_z;
        q1=altro_q1;
        q2=altro_q2;
        f_nearness=s_nearness;
    elseif (s_nearness>=f_nearness)
        allow=0;
    else
        z=altro_z;
        q1=altro_q1;
        q2=altro_q2;
        f_nearness=s_nearness;
    end
end
res=f_nearness;


% -----------------------------------------------
% function c_j_lu
% -----------------------------------------------
%
% Fast pseudo-solution of a least squares problem whose
% associated matrix is the Jacobian matrix of the
% gcd system (based on fast LU factorization).
%
% This function is used in the fast iterative refinement.
%


function [z]=c_j_lu(p,q,g,x,K)

% make sure that p,q,g are column vectors
[s1,s2]=size(p);
if s2>s1
    p=p.';
end
[s1,s2]=size(q);
if s2>s1
    q=q.';
end
[s1,s2]=size(g);
if s2>s1
    g=g.';
end

p=flipud(p);
q=flipud(q);
g=flipud(g);

% choose threshold used when checking whether
% the nullity of the Jacobian is too large

threshold=1e-9;

% compute degrees and other useful parameters
k=length(g)-1;
n=length(p)-1;
m=length(q)-1;
N=n+m+2*k+2+1;
M=n+m+k+3;
MN=lcm(N,M);
theta=exp(M*pi*i/MN);

% compute displacement generators for the Jacobian
zr=zeros(k,1);
G1=zeros(N,5);
G1(1,:)=[1 0 0 0 0];
G1(N,:)=[0 0 0 0 1];
G1(2:N-1,2)=[zr;p;zr;q(1:m)]-[g(2:k+1);zeros(N-k-2,1)];
G1(2:N-1,3)=[zeros(n,1);g;zeros(m+k,1)]-[zeros(n+k,1);g;zeros(m,1)];
G1(2:N-1,4)=[zeros(N-k-2,1);g(1:k)]-[theta*p(2:n+1); zr; theta*q;zeros(k,1,1)];
G2=zeros(5,M);
G2(2:4,[k+1, k+n+2, M])=eye(3);
G2(1,:)=[g(1:k).',g(k+1)-g(1),-p.',-q(1:m).',-q(m+1)-theta*p(1)];
G2(5,:)=[-g(2:k+1).', p(1)+q(m+1),p(2:n+1).', q.',g(k+1)-theta*g(1)];

% compute D0 (which is MxM), an auxiliary matrix
d0=ones(1,M);
a=0+(pi*i)/(MN);
alpha=exp(a);
for j=2:M
    d0(j)=d0(j-1)*alpha;
end
D0=diag(d0);

% compute D1 (which is NxN)
d1=ones(1,N);
beta=exp((2*pi*i)/N);
for j=2:N
    d1(j)=d1(j-1)*beta;
end
D1=d1;

% compute D2 (which is MxM and contains the
% shifted M-th roots of 1)
d2=ones(1,M);
beta=exp((2*pi*i)/M);
for j=2:M
    d2(j)=d2(j-1)*beta;
end
D2=exp(pi*i/MN)*d2;

% now D1 and D2 define a displacement operator for
% the Cauchy-like matrix C associated with the Jacobian

% compute displacement generators for this Cauchy-like matrix:
sN=sqrt(N);
sM=sqrt(M);
FD0=sM*ifft(D0);
G1_c=sN*ifft(G1);
G2_c=G2*FD0';

% now G1_c and G2_c are generators for C
% (i.e. D1*C - C*D2' = G1_c*G2_c)

% perform fast LU

L=zeros(N);
U=zeros(N,M);
P=eye(N);
P2=eye(M);


for k=1:M
    % the first generator is made orthogonal

    [Q,R]=qr(G1_c,0);
    G1_c=Q;
    G2_c=R*G2_c;

    % look for column of max magnitude in G2_c
    c=zeros(1,M-k+1);
    gg=abs(G2_c(:,k:M));
    gg=gg.*gg;
    c(k:M)=sum( gg,1);

    [val,col]=max(c);

    U(k,k:M)=(G1_c(k,:)*G2_c(:,k:M))./(D1(k)*ones(1,M-k+1)-D2(k:M));

    % swap s_k and s_col
    u=D2(k);
    D2(k)=D2(col);
    D2(col)=u;

    % swap psi_k and psi_col
    u=G2_c(:,k);
    G2_c(:,k)=G2_c(:,col);
    G2_c(:,col)=u;

    % swap columns in U
    u=U(:,k);
    U(:,k)=U(:,col);
    U(:,col)=u;

    % swap columns in P2
    u=P2(:,k);
    P2(:,k)=P2(:,col);
    P2(:,col)=u;
    L(k:N,k)=(G1_c(k:N,:)*G2_c(:,k))./((D1(k:N)-D2(k)*ones(1,N-k+1)).');
    [u,s]=max(abs(L(k:N,k)));
    s=s+k-1;
    U(k,k)=L(s,k);
    % swap t_k and t_q
    u=D1(k);
    D1(k)=D1(s);
    D1(s)=u;
    % swap phi_k and phi_q
    u=G1_c(k,:);
    G1_c(k,:)=G1_c(s,:);
    G1_c(s,:)=u;
    % swap rows in L
    u=L(k,:);
    L(k,:)=L(s,:);
    L(s,:)=u;
    % swap rows in P
    u=P(k,:);
    P(k,:)=P(s,:);
    P(s,:)=u;
    U(k,k+1:M)=(G1_c(k,:)*G2_c(:,k+1:M))./(D1(k)*ones(1,M-k)-D2(k+1:M));
    L(k,k)=1;
    L(k+1:N,k)=L(k+1:N,k)/U(k,k);
    G1_c(k+1:N,1)=G1_c(k+1:N,1)-G1_c(k,1)*L(k+1:N,k);
    G1_c(k+1:N,2)=G1_c(k+1:N,2)-G1_c(k,2)*L(k+1:N,k);
    G1_c(k+1:N,3)=G1_c(k+1:N,3)-G1_c(k,3)*L(k+1:N,k);
    G1_c(k+1:N,4)=G1_c(k+1:N,4)-G1_c(k,4)*L(k+1:N,k);
    G1_c(k+1:N,5)=G1_c(k+1:N,5)-G1_c(k,5)*L(k+1:N,k);
    G2_c(:,k+1:M)=G2_c(:,k+1:M)-G2_c(:,k)*(U(k,k+1:M)/U(k,k));

end
L(M+1:N,M+1:N)=eye(N-M);
%y=L\(P*sN*ifft([x;abs(norm(p)^2+norm(q)^2+norm(g)^2-K)]));
y=L\(P*sN*ifft([x;0]));

% detect nullity of Jacobian
nullity=0;
while abs(U(M-nullity,M-nullity))<threshold
    nullity=nullity+1;
end

y=y(1:M-nullity);
U=U(1:M-nullity,1:M-nullity);
y=U\y;
y=[y;zeros(nullity,1)];
z=D0'*(1/sM)*fft(P2*y);

% compute a vector in Ker(S)
nsp=[g; -p; -q];
nsp=nsp/norm(nsp);
z=z-dot(nsp,z)*nsp;   % projection of z on range(S)


% --------------------------------------------------
% function c_cofactors
% --------------------------------------------------
%
% Uses a fast algorithm (GKO + almost complete pivoting)
% to compute the cofactors of polynomials f and g
% with respect to an approximate gcd of degree tdeg:
%
% [f1,f2]=c_cofactors(f,g,tdeg);
%

function [f1,f2]=c_cofactors(f,g,tdeg)

threshold=1e-9;  % tolerance for rank determination

% make sure that f and g are column vectors
[s1,s2]=size(f);
if s2>s1
    f=f.';
end
[s1,s2]=size(g);
if s2>s1
    g=g.';
end

% reverse polynomial for agreement with other programs
f=flipud(f);
g=flipud(g);

% compute degrees and other useful parameters
n=length(f)-1;
m=length(g)-1;
N=n+m-tdeg+1;
M=n+m-2*tdeg+2;
MN=lcm(N,M);
theta=exp(M*pi*i/MN);

% compute displacement generators for the cofactor matrix
G1=zeros(N,2);
G1(:,1)=[f(n+1); zeros(m-tdeg,1); f(1:n)]-[g; zeros(n-tdeg,1)];
G1(:,2)=[g(m+1); zeros(n-tdeg,1); g(1:m)]-theta*[f; zeros(m-tdeg,1)];
G2=zeros(2,M);
G2(:,[m-tdeg+1,m+n-2*tdeg+2])=eye(2);

% compute D0 (which is MxM), an auxiliary matrix
d0=ones(1,M);
a=0+(pi*i)/(MN);
alpha=exp(a);
for j=2:M
    d0(j)=d0(j-1)*alpha;
end
D0=diag(d0);

% compute D1 (which is NxN)
d1=ones(1,N);
beta=exp((2*pi*i)/N);
for j=2:N
    d1(j)=d1(j-1)*beta;
end
D1=d1;

% compute D2 (which is MxM and contains the shifted M-th
% roots of 1)
d2=ones(1,M);
beta=exp((2*pi*i)/M);
for j=2:M
    d2(j)=d2(j-1)*beta;
end
D2=exp(pi*i/MN)*d2;

% now D1 and D2 define a displacement operator for the
% Cauchy-like matrix C associated with the cofator matrix

% compute displacement generators G1_c and G2_c for C
sN=sqrt(N);
FD0=sN*ifft(D0);
G1_c=sN*ifft(G1);
G2_c=G2*FD0';

% perform fast LU
L=zeros(N);
U=zeros(N,M);

P2=eye(M);

for k=1:M
    % the first generator is made orthogonal
    [Q,R]=qr(G1_c,0);
    G1_c=Q;
    G2_c=R*G2_c;

    % look for column of max magnitude in G2_c
    c=zeros(1,M-k+1);
    c(k:M)=conj(G2_c(1,k:M)).*G2_c(1,k:M)+conj(G2_c(2,k:M)).*G2_c(2,k:M);

    [val,col]=max(c);

    U(k,k:M)=(G1_c(k,:)*G2_c(:,k:M))./((D1(k)*ones(1,M-k+1)-D2(k:M)));


    % swap s_k and s_col
    u=D2(k);
    D2(k)=D2(col);
    D2(col)=u;

    % swap psi_k and psi_col
    u=G2_c(:,k);
    G2_c(:,k)=G2_c(:,col);
    G2_c(:,col)=u;

    % swap columns in U
    u=U(:,k);
    U(:,k)=U(:,col);
    U(:,col)=u;

    % swap columns in P2
    u=P2(:,k);
    P2(:,k)=P2(:,col);
    P2(:,col)=u;
    L(k:N,k)=(G1_c(k:N,:)*G2_c(:,k))./((D1(k:N)-D2(k)*ones(1,N-k+1)).');


    [u,s]=max(abs(L(k:N,k)));
    s=s+k-1;
    U(k,k)=L(s,k);
    % swap t_k and t_q
    u=D1(k);
    D1(k)=D1(s);
    D1(s)=u;
    % swap phi_k and phi_q
    u=G1_c(k,:);
    G1_c(k,:)=G1_c(s,:);
    G1_c(s,:)=u;
    % swap rows in L
    u=L(k,:);
    L(k,:)=L(s,:);
    L(s,:)=u;


    U(k,k+1:M)=(G1_c(k,:)*G2_c(:,k+1:M))./(D1(k)*ones(1,M-k)-D2(k+1:M));


    L(k,k)=1;
    L(k+1:N,k)=L(k+1:N,k)/U(k,k);
    G1_c(k+1:N,1)=G1_c(k+1:N,1)-G1_c(k,1)*L(k+1:N,k);
    G1_c(k+1:N,2)=G1_c(k+1:N,2)-G1_c(k,2)*L(k+1:N,k);

    G2_c(:,k+1:M)=G2_c(:,k+1:M)-G2_c(:,k)*(U(k,k+1:M)/U(k,k));

end

% compute a vector in the null space of the cofactor matrix
emme=m+n-2*tdeg+2;
nullity=0;
while abs(U(emme-nullity,emme-nullity))<threshold
    nullity=nullity+1;
end
if nullity==0
    nullity=1;
end
nullity=nullity-1; % this is not the nullity now
r=-U(1:emme-1-nullity,emme-nullity);
U=U(1:emme-1-nullity,1:emme-1-nullity);
x=U\r;
x=[x;1;zeros(nullity,1)];
x=D0'*fft(P2*x);  % this vector contains the cofactors
f1=flipud((x(1:m-tdeg+1)));
f2=flipud((x(m-tdeg+2:emme)));


% -------------------------------------------------
% function c_gs_gecp
% -------------------------------------------------
%
% Fast LU factorization for a Sylvester matrix.
%
% The function only returns the diagonal (in
% absolute value) of the upper triangular
% factor U, since it contains all the required information.
%

function dU = c_gs_gecp(f,g)

% compute degrees of input polynomials
n=length(f)-1;
m=length(g)-1;
N=m+n;

% compute displacement generators for the
% Sylvester matrix
G2=zeros(2,N);
G2(1,N-m:N)=g;
G2(1,1:n)=G2(1,1:n)-f(2:n+1);
G2(1,N)=G2(1,N)+f(1);
G2(2,N-n:N)=f;
G2(2,1:m)=G2(2,1:m)-g(2:m+1);
G2(2,N)=G2(2,N)+g(1);
G1=zeros(N,2);
G1([1;m+1],:)=eye(2);

% compute the auxiliary matrix D0
d0=ones(1,N);
a=0+(pi*i)/N;
alpha=exp(a);
for j=2:N
    d0(j)=d0(j-1)*alpha;
end
D0=diag(d0);

% compute D1
d1=ones(1,N);
beta=alpha*alpha;
for j=2:N
    d1(j)=d1(j-1)*beta;
end
D1=d1;
% compute D2
D2=alpha*d1;
% now D1 and D2 are the diagonals of the
% matrices that define a displacement operator for
% the Cauchy-like matrix C associated with
% the Sylvester matrix

% compute generators G1_c and G2_c  for C
sN=sqrt(N);
FD0=sN*ifft(D0);

G1_c=sN*ifft(G1);
G2_c=G2*FD0';

% perform fast LU

L=zeros(N);
U=zeros(N);

for k=1:N
    % the first generator is made orthogonal
    [Q,R]=qr(G1_c,0);
    G1_c=Q;
    G2_c=R*G2_c;

    % look for column of max magnitude in G2_c
    c=zeros(1,N-k+1);
    c(k:N)=conj(G2_c(1,k:N)).*G2_c(1,k:N)+conj(G2_c(2,k:N)).*G2_c(2,k:N);

    [val,col]=max(c);
    U(k,k:N)=(G1_c(k,:)*G2_c(:,k:N))./(D1(k)-D2(k:N));

    % swap s_k and s_col
    u=D2(k);
    D2(k)=D2(col);
    D2(col)=u;

    % swap psi_k and psi_col
    u=G2_c(:,k);
    G2_c(:,k)=G2_c(:,col);
    G2_c(:,col)=u;

    % swap columns in U
    u=U(:,k);
    U(:,k)=U(:,col);
    U(:,col)=u;

    L(k:N,k)=(G1_c(k:N,:)*G2_c(:,k))./((D1(k:N)-D2(k)*ones(1,N-k+1)).');


    [u,s]=max(abs(L(k:N,k)));
    s=s+k-1;
    U(k,k)=L(s,k);
    % swap t_k and t_q
    u=D1(k);
    D1(k)=D1(s);
    D1(s)=u;
    % swap phi_k and phi_q
    u=G1_c(k,:);
    G1_c(k,:)=G1_c(s,:);
    G1_c(s,:)=u;
    % swap rows in L
    u=L(k,:);
    L(k,:)=L(s,:);
    L(s,:)=u;


    U(k,k+1:N)=(G1_c(k,:)*G2_c(:,k+1:N))./(D1(k)*ones(1,N-k)-D2(k+1:N));


    L(k,k)=1;

    L(k+1:N,k)=L(k+1:N,k)./U(k,k);
    G2_c(:,k+1:N)=G2_c(:,k+1:N)-G2_c(:,k)*(U(k,k+1:N)./U(k,k));
    G1_c(k+1:N,1)=G1_c(k+1:N,1)-G1_c(k,1)*L(k+1:N,k);
    G1_c(k+1:N,2)=G1_c(k+1:N,2)-G1_c(k,2)*L(k+1:N,k);


end
dU=abs(diag(U));

% -----------------------------------------------------
% function c_iterfast
% -----------------------------------------------------
%
% writes the system associated with the GCD
%

function newz=c_iterfast(a,b,q1,q2,g,K,optflag)

threshold=1e-6;
c=flipud([conv(q2,g).';conv(q1,g).']);
direction=c_j_lu(q1,q2,g,(c-[(fliplr(a)).'; (fliplr(b)).']),K);
c1=costfun(1,a,b,q1,q2,g,direction);
if (optflag==1)&&(c1>costfun(0,a,b,q1,q2,g,direction))&&(c1>threshold)
    %if (optflag==1)&&(c1>costfun(0,a,b,q1,q2,g,direction))
    options=optimset('Display','off','MaxIter',3);
    fprintf('*\n')
    alpha=fminbnd(@(apar) costfun(apar,a,b,q1,q2,g,direction),0,5,options);
else
    alpha=1;
end
newz=[(fliplr(g).'); fliplr(q1).'; fliplr(q2).'] - alpha*direction;


% -----------------------------------------------------
% function costfun
% -----------------------------------------------------
%
% computes the cost function (to be used for line search
% in c_iterfast)
%

function F=costfun(apar,a,b,q1,q2,g,direction)

z_t=[(fliplr(g).'); fliplr(q1).'; fliplr(q2).'] - apar*direction;
lg=length(g);
lq1=length(q1)+lg;
newg=flipud(z_t(1:lg));
polya=conv(flipud(z_t(lg+1:lq1)),newg);
polyb=conv(flipud(z_t(lq1+1:end)),newg);
sz=size(polya);
if sz(1)>sz(2)
    polya=polya.';
    polyb=polyb.';
end
F=norm([a-polya, b-polyb]);


% --------------------------------------------------------
% function div_fft
% --------------------------------------------------------
%
% Computes polynomial division using the FFT
% (that is, with an evaluation/interpolation technique)
%
% u=div_fft(f,g)
% where u=f/g.
%
% Notice that:
% 1) f must be exactly divisible by g;
% 2) if a root of 1 ia also a root of g, the function
%    attempts scaling of f and g and division with
%    evaluation/iterpolation;
% 3) if even this attmpt fails, the function
%    uses lsdiv instead of the FFT.
%
% Last modified 27 October 2008.
%

function u=div_fft(f,g,delta)
if delta>5e-7
    u=c_lsdiv(f,g);
else
    tol=1e-8;
    n=length(f);
    m=length(g);
    % controllo che f e g siano vettori riga
    s=size(f);
    if s(1)>s(2)
        f=f.';
    end
    s=size(g);
    if s(1)>s(2)
        g=g.';
    end

    k=n-m+1; % length(u)
    c=ceil(n/k);
    r=k*c;
    fpad=[zeros(1,r-n) f];
    gpad=[zeros(1,r-m) g];
    fval=fft(fpad);
    gval=fft(gpad);
    if min(abs(gval))>tol
        uval=fval(1:c:r)./gval(1:c:r);
        uu=ifft(uval);
        u=zeros(1,k);
        u(k)=uu(1);
        u(1:k-1)=uu(2:k);
    else
        alpha=1+1/n;
        scale=ones(1,n);
        for j=1:n-1
            scale(j+1)=alpha*scale(j);
        end
        fs=f.*scale;
        gs=g.*scale(1:m);
        fpad=[zeros(1,r-n) fs];
        gpad=[zeros(1,r-m) gs];
        fval=fft(fpad);
        gval=fft(gpad);
        if min(abs(gval))>tol
            uval=fval(1:c:r)./gval(1:c:r);
            uu=ifft(uval);
            u=zeros(1,k);
            u(k)=uu(1);
            u(1:k-1)=uu(2:k);
            u=u./scale(1:k);
        else
            fprintf('Using c_lsdiv...\n')
            u=c_lsdiv(f,g);
        end
    end
end

%-------------------------------------------------
% function c_lsdiv
%-------------------------------------------------
%
% Computes polynomial division by solving a
% linear least squares system with GKO.
% The function computes u=f/g.
%
% (Bug fixed on October 9, 2010; thanks to 
% Jean-Philippe Chancelier for pointing it out).

function u=c_lsdiv(f,g)
sf=size(f);
if sf(1)<sf(2)
    f=f.';
    g=g.';
end
n=length(f)-1;
m=length(g)-1;
k=n-m;
threshold=1e-13;

N=n+1;
M=k+1;  % size of convolution matrix
MN=lcm(N,M);
theta=exp(M*pi*i/MN);

% choose Toeplitz-like generators
G1=zeros(n+1,1);
G1(1)=g(m+1);
G1(k+2:n+1)=g(1:m);
G1(1:m+1)=G1(1:m+1)-theta*g;
G2=[zeros(1,k),1];

d0=ones(1,M);
a=0+(pi*i)/(MN);
alpha=exp(a);
for j=2:M
    d0(j)=d0(j-1)*alpha;
end
D0=diag(d0);

% compute D1 (which is NxN)
d1=ones(1,N);
beta=exp((2*pi*i)/N);
for j=2:N
    d1(j)=d1(j-1)*beta;
end
D1=d1;

% compute D2 (which is MxM and contains the shifted M-th
% roots of 1)
d2=ones(1,M);
beta=exp((2*pi*i)/M);
for j=2:M
    d2(j)=d2(j-1)*beta;
end
D2=exp(pi*i/MN)*d2;

% now D1 and D2 define a displacement operator for the
% Cauchy-like matrix C associated with the cofator matrix

% compute displacement generators G1_c and G2_c for C
sN=sqrt(N);
sM=sqrt(M);
FD0=sN*ifft(D0);
G1_c=sN*ifft(G1);
G2_c=G2*FD0';

% perform fast LU
L=zeros(N);
U=zeros(N,M);
P=eye(N);
P2=eye(M);

for k=1:M
    % the first generator is made orthogonal
    [Q,R]=qr(G1_c,0);
    G1_c=Q;
    G2_c=R*G2_c;

    % look for column of max magnitude in G2_c
    c=zeros(1,M-k+1);
    c(k:M)=conj(G2_c(k:M)).*G2_c(k:M);

    [val,col]=max(c);

    U(k,k:M)=(G1_c(k)*G2_c(k:M))./((D1(k)*ones(1,M-k+1)-D2(k:M)));


    % swap s_k and s_col
    u=D2(k);
    D2(k)=D2(col);
    D2(col)=u;

    % swap psi_k and psi_col
    u=G2_c(k);
    G2_c(k)=G2_c(col);
    G2_c(col)=u;

    % swap columns in U
    u=U(:,k);
    U(:,k)=U(:,col);
    U(:,col)=u;

    % swap columns in P2
    u=P2(:,k);
    P2(:,k)=P2(:,col);
    P2(:,col)=u;
    L(k:N,k)=(G1_c(k:N)*G2_c(k))./((D1(k:N)-D2(k)*ones(1,N-k+1)).');


    [u,s]=max(abs(L(k:N,k)));
    s=s+k-1;
    U(k,k)=L(s,k);
    % swap t_k and t_q
    u=D1(k);
    D1(k)=D1(s);
    D1(s)=u;
    % swap phi_k and phi_q
    u=G1_c(k);
    G1_c(k)=G1_c(s);
    G1_c(s)=u;
    % swap rows in L
    u=L(k,:);
    L(k,:)=L(s,:);
    L(s,:)=u;
    % swap rows in P
    u=P(k,:);
    P(k,:)=P(s,:);
    P(s,:)=u;


    U(k,k+1:M)=(G1_c(k)*G2_c(k+1:M))./(D1(k)*ones(1,M-k)-D2(k+1:M));


    L(k,k)=1;
    L(k+1:N,k)=L(k+1:N,k)/U(k,k);
    G1_c(k+1:N)=G1_c(k+1:N)-G1_c(k)*L(k+1:N,k);

    G2_c(k+1:M)=G2_c(k+1:M)-G2_c(k)*(U(k,k+1:M)/U(k,k));

end
L(M+1:N,M+1:N)=eye(N-M);
y=L\(P*sN*ifft(f));

% detect nullity
nullity=0;
while abs(U(M-nullity,M-nullity))<threshold
    nullity=nullity+1;
end

y=y(1:M-nullity);
U=U(1:M-nullity,1:M-nullity);
y=U\y;
y=[y;zeros(nullity,1)];
u=D0'*(1/sM)*fft(P2*y);
%if sf(1)<sf(2)
u=u.';
%end