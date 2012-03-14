
# Stubs to call STLS
slra <- function(P, S, R = (total_cols - 1), X0 = NULL, 
            opt = list(), compute.dp = FALSE, ret.obj=FALSE) {

  ###### Parse structure
  # convert stucture into list
  if (!is.list(S)) {
    S <- list(1,S);
    names(S) <- c("k", "A");
  } 
  storage.mode(S$k) <- 'integer';
  if (S$k <= 0) {
    stop ('Incorrect row dimension of the block');
  }
  if (is.vector(S$A)) {
    S$A <- matrix(S$A, 1, length(S$A));
  }

  # Parse structure matrix A
  storage.mode(S$A) <- 'double';
  if (nrow(S$A) > 10) {
    stop("There is more than 10 blocks");
  }
  if (ncol(S$A) > 4) {
     stop("Structure matrix has incorrect number of columns");
  }
  if (ncol(S$A) < 2) {
    S$A <- cbind(S$A, matrix(1, nrow(S$A), 1));
  }

  total_cols <- sum(S$A[,1] * S$A[,2]);
  storage.mode(total_cols) <- 'integer';

  if (is.vector(X0)) {
    X0 <- matrix(X0, length(X0), 1);
  }
   
  storage.mode(R) <- 'integer';
  n <- R;
  d <- total_cols - n;   

  if (!is.null(X0)) {
    storage.mode(X0) <- 'double';
  
    if((nrow(X0) != n || ncol(X0) != d)) {
      stop('nrow(X0) != n || ncol(X0) != d');
    }
  }

  sobj <- .Call("create_slra_object", n, d, P, S, X0, opt, compute.dp);
  
  if (ret.obj) {
    res <- sobj; 
  } else {
    .Call("optimize_slra", sobj);
    res <- .Call("return_slra_object", sobj);
    if (compute.dp) {
      res <- c(res, list(ph=.Call("compute_slra_correction", sobj)));
    }
  }

  res;    
}

optimize_slra <- function(sobj) {
  .Call("optimize_slra", sobj);
}

evaluate_slra_f <- function(sobj, ...) {
  .Call("evaluate_slra_f", sobj, ...);
}

evaluate_slra_df <- function(sobj, ...) {
  .Call("evaluate_slra_df", sobj, ...);
}

return_slra_object <- function(sobj) {
  .Call("return_slra_object", sobj);
}

compute_slra_correction <- function(sobj) {
  .Call("compute_slra_correction", sobj);
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




