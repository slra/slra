
# Stubs to call STLS
slra <- function(P, S, R = (total_cols - 1), X = NULL, 
            opts = list(epsabs = 0, epsrel = 1e-5, epsgrad = 1e-5, 
                        maxiter = 100, disp = 'iter'),
            compute.dp = FALSE) {

  ###### Parse structure
  # convert stucture into list
  if (!is.list(S)) {
    S <- list(1,S);
    names(S) <- c("k", "A");
  } 
  storage.mode(S$k) <- 'integer'
  if (S$k <= 0) {
    stop ('Incorrect row dimension of the block');
  }
  if (is.vector(S$A)) {
    S$A <- matrix(S$A, 1, length(S$A));
  }

  # Parse structure matrix A
  storage.mode(S$A) <- 'integer'
  if (nrow(S$A) > 10) {
    stop("There is more than 10 blocks");
  }
  if (ncol(S$A) > 4) {
     stop("Structure matrix has incorrect number of columns");
  }
  if (!prod(S$A[,1] >= 0)) {
     stop("Structure matrix contains negative numbers");
  }
  if (ncol(S$A) < 2) {
    S$A <- cbind(S$A, matrix(1, nrow(S$A), 1));
  }
  if (ncol(S$A) < 4) {
    S$A <- cbind(S$A, matrix(0, nrow(S$A), 4 - ncol(S$A)));
  }
  





  if (!prod((S$A[,1] >= 1) & (S$A[,1] <= 4))) {
    stop('Unrecognized block structure');
  }

  total_cols <- sum(S$A[,2]);

  if (is.vector(X)) {
    X <- matrix(X, length(X), 1);
  }
   
  if (is.null(P)) {
    stop("P should be given");
  }

  storage.mode(R) <- 'integer';
  n <- R;
  d <- total_cols - n;   

  if (!is.null(X)) {
    storage.mode(X) <- 'double';
  
    if((nrow(X) != n || ncol(X) != d)) {
      stop('nrow(X) != n || ncol(X) != d');
    }
  }


  .Call("rslra", n, d, P, S, X, opts, compute.dp);
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




