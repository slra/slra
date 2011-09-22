% constructs the rhs for the KKT equations that arise in the STLN 
% optimization problem for Toeplitz matrices [A b]

function [rhs,p2] = constrhs(alpha,beta,x,lambda,t,m,n)

p11 = [-alpha;-beta;zeros(n,1)];   

%construction of p12
perm1 = [2:m+n 1 m+n+1:m+2*n]'; %permutation for ease of computations
xspec = [-x.' 1];
   jlp = zeros(m+2*n,1);
for i = 1:n
   jlp(i,1) = xspec(n+2-i:n+1)*lambda(1:i);
end
for i = n+1:m
   jlp(i,1) = xspec*lambda(i-n:i);
end

for i = m+1:m+n
   jlp(i,1) = xspec(1:m+n-i+1)*lambda(i-n:m);
end

jlp(m+n+1:m+2*n,1) = -toeplitz(t(n+1:n+m)+alpha(n:n+m-1),...
    t(n+1:-1:2)+alpha(n:-1:1)).'*lambda;
p12 = jlp(perm1);

%construction of p2
rt = t(1:m)-toeplitz(t(n+1:n+m)+alpha(n:n+m-1),t(n+1:-1:2)+alpha(n:-1:1))*x;
rt(1) = rt(1)+beta;
rt(2:m) = rt(2:m)+alpha(1:m-1);
p2 = rt;

rhs = [p11+p12;p2];


