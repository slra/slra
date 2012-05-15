library('pracma')
library('DAAG')
library('Rslra');

tls <- function(M, d) {
  nd <- dim(M)[2];
  R_tls <- t(svd(M)$v[,(nd-d+1):nd]);
  if (d == 1) {
    res <- -R_tls[1:(nd-d)] / R_tls[nd];
  } else {
    res <- (-qr.solve(R_tls[1:d,(nd-d+1):nd], R_tls[1:d,1:(nd-d)]));
  }
  res;
}


cat("\n0) Hankel low-rank approximation example:\n");
N <- 33
K <- 3
L <- N-K+1
f0 <- sin(1:N * pi /12 )
fn <- rnorm(N, 0, 0.2)
f <- f0 + fn

h0 <- outer(1:L, 1:K, function(x,y) f0[x+y-1]);
h <- outer(1:L, 1:K, function(x,y) f[x+y-1]);

amat <- h[,1:(K-1)]
bmat <- h[,K]
  
res <- slra(f, list(m=K), K-1, opt = list(disp='iter'), compute.Rh=TRUE)
cat('X0 = ', tls(h0,1), "\n");
cat('Xtls = ', tls(h,1), "\n");
cat('Xslra = ', res$info$Rh[1:(K-1)], "\n\n");
pause();

# ===================================================
# ======                                        =====
# ======              STLS demo                 =====
# ======                                        =====
# ===================================================
# ===================================================
#

# ===================================================
cat("\n1) First we solve a least squares problem\n")

# Define dimensions and generate random data
m <- 100; n <- 5; d <- 2;
a <- matrix(runif(m*n),m,n);
b <- matrix(runif(m*d),m,d);

# Find the LS estimate by Matlab's \
system.time(x_ls <- qr.solve(a, b))

# Define and solve the LS problem as a (very special) STLS problem
s_ls <- list(m = rep(1, n+d), n = m, w = c(rep(Inf, n), rep(1, d)));
print(system.time(res_ls <- slra(c(as.vector(a),as.vector(b)), s_ls, n,
                                 opt = list(disp='iter'))))
cat("x_ls = \n");
print(t(x_ls));
cat("x_slra = \n");
print(res_ls$info$Rh[1:d,1:n]);
cat("\n")
pause();


# ===================================================
cat("\n2) Next we solve a TLS problem with (TLS+LS)/2 initial approximation\n")

# Define and solve the TLS problem as an STLS problem
s_tls <- list(m = rep(1, n+d), n = m);
X_tls <- tls(cbind(a,b), d);
R_tls <- cbind(X_tls,-eye(d));
print(system.time(res <- slra(c(as.vector(a),as.vector(b)),s_tls,n, 
                            opt = list(Rini = (res_ls$info$Rh +R_tls)/2,
                            maxiter=1000, disp='iter'))))
cat("X_tls = \n")
print(R_tls[1:d,1:n]);
cat("X_slra = \n")
print(res$info$Rh[1:d,1:n]);
cat("\n");
# Press any key to continue
pause();


cat("\n2) Next we solve a deconvolution problem\n")
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
v_n  <- 0.05; # noise level
p_a  <- p_a0 + v_n * rnorm(n+m-1);
a    <- outer(1:m, 1:n, function(x,y) p_a[x-y+n]);
b    <- b0 + v_n * rnorm(m);

# Ignore the structure and estimate via LS and TLS
#print(
system.time(xh_ls <- qr.solve(a,b))
#);
#print(system.time(xh_tls <- tls(a,b)));

# Define the structure and solve the deconvolution problem via STLS
s <- list(m=c(2,1), phi=matrix(c(0, 1, 0,
                                 1, 0, 0,
                                 0, 0, 1), 3, 3));
#print(
system.time(res <- slra(c(p_a,b), s, n))
xh_tls <- tls(cbind(a,b),1)
xh_slra <- res$info$Rh[1:n];
#); 
cat("X_ls = ", xh_ls, "\n");
cat("X_ls = ", xh_tls, "\n");
cat("X_slra = ", xh_slra, "\n");
#print(res$info$fmin) # value of the cost function at xh_stls

# Solve via an alternative STLS method
#tic, xh_stln <- faststln1(a,b); t_stln <- toc
#xh_stln(1:2)'
#cost1(xh_stln,a,b,s) # value of the cost function at xh_stln

norm2 <- function(x) { sqrt(sum(x*x)); }

# Compare the relative errors of estimation 
cat("e_ls = ", e_ls   <- norm2(xh_ls-x0)/norm2(x0), "\n"); 
cat("e_tls = ",e_tls  <- norm2(xh_tls-x0)/norm2(x0), "\n");  
cat("e_slra = ", e_stls <- norm2(xh_slra-x0)/norm2(x0), "\n");   

# Press any key to continue
pause()





