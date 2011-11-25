
# Stubs to call STLS
stls <- function(A = NULL, B = NULL, S, X = NULL, opts = list(epsabs = 0, epsrel = 1e-6, epsgrad = 1e-6, maxiter = 100, disp = 'iter'),
           P = NULL, compute.dp = !is.null(P)) {

  ###### Parse structure
  if (!is.list(S)) {
    S <- list(1,S);
    names(S) <- c("k", "A");
  } 
  
  storage.mode(S$k) <- 'integer'
  if (!is.null(S$m)) {
    storage.mode(S$m) <- 'integer'
  }
  if (!is.null(S$d)) {
    storage.mode(S$d) <- 'integer'
  }

  if (S$k <= 0) {
    stop ('Incorrect row dimension of the block');
  }
  if (is.vector(S$A)) {
    S$A <- matrix(S$A, 1, length(S$A));
  }
  
  storage.mode(S$A) <- 'integer'
  if (nrow(S$A) > 10) {
    stop("There is more than 10 blocks");
  }
  if (ncol(S$A) != 3) {
     stop("Structure specification matrix has incorrect number of columns (!= 3)");
  }
  if (!prod((S$A[,1] >= 1) & (S$A[,1] <= 4))) {
    stop('Unrecognized block structure');
  }

  n_plus_d <- sum(S$A[,2]);
    

  if (is.vector(X)) {
    X <- matrix(X, length(X), 1);
  }

   
  if (is.null(A) || is.null(B)) {
    if (is.null(P)) {
      stop("At least one of (A,B) or P should be given");
    }
  
    if (!is.matrix(X)) {
      if (is.null(S$d)) {
        stop("d is not defined");
      }
      d <- S$d;
      n <- n_plus_d - d;
    } else {
      n <- nrow(X);
      d <- ncol(X);
    }
  } else {
    if (is.vector(A)) {
      A <- matrix(A, length(A),1);
    } 
    if (is.vector(B)) {
      B <- matrix(B, length(B),1);
    } 

    m <- nrow(A);
    n <- ncol(A);
    d <- ncol(B);

    storage.mode(A) <- storage.mode(B) <- 'double';
    if (m != nrow(B)) {
      stop('m != nrow(B)');
    }
  }

  if (!is.null(X)) {
    storage.mode(X) <- 'double';
  
    if((nrow(X) != n || ncol(X) != d)) {
      stop('nrow(X) != n || ncol(X) != d');
    }
  }


  .Call("rstls", n, d, A, B, S, X, opts, P, compute.dp);
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




