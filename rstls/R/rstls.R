
# Stubs to call STLS
stls <- function(A, B, S, X=tls(A,B), opts = list(epsabs = 0, epsrel = 1e-6, epsgrad = 1e-6, maxiter = 100, disp = 'iter')) {
  if (!is.matrix(A)) {
    stop('A is not a matrix');
  }

  if (!is.matrix(B)) {
     B <- matrix(B, length(B),1);
  } 

  if (!is.matrix(X)) {
     X <- matrix(X, length(X), 1);
  } 

  m <- nrow(A);
  n <- ncol(A);
  d <- ncol(B);

  if (m != nrow(B)) {
     stop('m != nrow(B)');
  }

  if (nrow(X) != n || ncol(X) != d) {
     stop('nrow(X) != n || ncol(X) != d');
  }
  

  storage.mode(S$k) <- 'integer'
  if (S$k <= 0 || m %% S$k != 0) {
    stop ('Incorrect row dimension of the block');
  }

  storage.mode(S$A) <- 'integer'
  if (nrow(S$A) > 10 || ncol(S$A) != 3) {
    stop ('Incorrect structure matrix');
  }

  if (!prod((S$A[,1] >= 1) & (S$A[,1] <= 4))) {
    stop('Unrecognized block structure');
  }

  storage.mode(A) <- storage.mode(B) <- storage.mode(X) <- 'double'

  .Call("rstls", A, B, S, X, opts);
}


tls <- function(A, B) {
  if (!is.matrix(A)) {
    stop('A is not a matrix');
  }

  if (!is.matrix(B)) {
     B <- matrix(B, length(B),1);
  } 

  if (nrow(A) != nrow(B)) {
     stop('nrow(A) != nrow(B)');
  }


  storage.mode(A) <- storage.mode(B) <- 'double'

  .Call("rtls", A, B);
}




