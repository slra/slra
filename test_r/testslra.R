############################
# Hankel lra example
############################
library('DAAG')
library('Rslra');

N <- 33
K <- 3
L <- N-K+1
f0 <- sin(1:N * pi /12 )
fn <- rnorm(N, 0, 0.2)
f <- f0 + fn

h0 <- outer(1:L, 1:K, function(x,y) f0[x+y-1]);
h <- outer(1:L, 1:K, function(x,y) f[x+y-1]);

xext0 <- svd(h0)$v[,3]
xext <- svd(h)$v[,3]
x0 <- -xext0[1:2]/xext0[3]
x <- -xext[1:2]/xext[3]

amat <- h[,1:2]
bmat <- h[,3]
xmat <- matrix(x, 2, 1)

res <- slra(f, matrix(c(2, 3, 1), 1, 3), 2, xmat, opt = list(method='l'))
print('x0(exact), x(tls), x(stls)')
print(x0)
print(x)
print(res)

pause();

# ===================================================
# ======                                        =====
# ======              STLS demo                 =====
# ======                                        =====
# ===================================================
# ===================================================
#
#
# ===================================================
# ====== First we solve a least squares problem =====
# ===================================================
#
# Define dimensions and generate random data
m <- 100; n <- 5; d <- 2;
a <- matrix(runif(m*n),m,n);
b <- matrix(runif(m*d),m,d);

# Find the LS estimate by Matlab's \
system.time(x_ls <- qr.solve(a, b))
print(x_ls)

# Define and solve the LS problem as a (very special) STLS problem
s_ls <- rbind(c(4, n, 1), c(3, d, 1));
pause();
print(system.time(res <- slra(c(as.vector(t(a)),as.vector(t(b))), s_ls, n)))
print(res$xh)

# Press any key to continue 
pause();


# ===================================================
# ======      Next we solve a TLS problem       =====
# ===================================================
#
# The data is a,b used above.
#
# Solve the TLS problem via SVD
system.time(x_tls <- tls(a,b));
print(x_tls);
pause();

# Define and solve the TLS problem as an STLS problem
s_tls <- list(k=1, A=matrix(c(3, n+d, 1), 1,3));
print(system.time(res <- slra(as.vector(t(cbind(a,b))),s_tls,n)))
print(res$xh)

# Press any key to continue
pause();


# ===================================================
# ======  Next we solve a deconvolution problem =====
# ===================================================
#
# b0 <- conv(p_a0,x0), p_a <- p_a0 + noise, b <- b0 + noise
# Problem: given a and b, estimate x0

m <- 200; # length(p_a0)
n <- 2;   # length(b0)

# Generate true data: p_a0 and b0
p_a0 <- runif(n+m-1);
a0   <- outer(1:m, 1:n, function(x,y) p_a0[x-y+n]);
x0   <- runif(n);
b0   <- a0 %*% x0;

# Add noise: p_a <- p_a0 + noise, b <- b0 + noise
v_n  <- 0.25; # noise level
p_a  <- p_a0 + v_n * rnorm(n+m-1);
a    <- outer(1:m, 1:n, function(x,y) p_a[x-y+n]);
b    <- b0 + v_n * rnorm(m);

# Ignore the structure and estimate via LS and TLS
print(system.time(xh_ls <- qr.solve(a,b)));
print(system.time(xh_tls <- tls(a,b)));

# Define the structure and solve the deconvolution problem via STLS
s <- list(k=1, A=rbind(c(1, n, 1), c(3, 1, 1)));
print(system.time(res <- slra(c(p_a,b), s, n))); 
print(xh_stls <- res$xh)
print(res$info$fmin) # value of the cost function at xh_stls

# Solve via an alternative STLS method
#tic, xh_stln <- faststln1(a,b); t_stln <- toc
#xh_stln(1:2)'
#cost1(xh_stln,a,b,s) # value of the cost function at xh_stln

norm2 <- function(x) { sqrt(sum(x*x)); }

# Compare the relative errors of estimation 
cat("e_ls = ", e_ls   <- norm2(xh_ls-x0)/norm2(x0), "\n"); 
cat("e_tls = ",e_tls  <- norm2(xh_tls-x0)/norm2(x0), "\n");  
cat("e_stls = ", e_stls <- norm2(xh_stls-x0)/norm2(x0), "\n");   

# Press any key to continue
pause()





