 
#  storage.mode(R) <- 'integer';
#  n <- R;
#  d <- total_cols - n;   

#  if (!is.null(X0)) {
#    storage.mode(X0) <- 'double';
#  
#    if((nrow(X0) != n || ncol(X0) != d)) {
#      stop('nrow(X0) != n || ncol(X0) != d');
#    }
#  }

#  sobj <- .Call("create_slra_object", n, d, P, S, X0, opt, compute.dp);
  
#  if (ret.obj) {
#    res <- sobj; 
#  } else {
#    .Call("optimize_gsl", sobj, opt);
#    res <- .Call("return_slra_object", sobj);
#    if (compute.dp) {
#      res <- c(res, list(ph=.Call("compute_slra_correction", sobj)));
#    }
#  }
#  res;    
#}


slra <- function(p, s, r = sum(S$m) - 1, 
            opt = list(), compute.dp = FALSE, ret.obj=FALSE) {
  # Check necessary parameters            
  if (!is.list(s) || is.null(s$m) || is.null(s$n)) {
    stop ('Structure must be a list with "m" and "n" elements');
  } 
  storage.mode(s$m) <- storage.mode(s$n) <- 'integer';
  storage.mode(s$m) <- storage.mode(s$n) <- 'double';
  if (!prod(s$m > 0) || !prod(s$n > 0)) {
    stop ('"m" and "n" elements should be ');
  }
  
  
  
  
  storage.mode(r) <- 'integer';
  storage.mode(compute.dp) <- 'integer';

  res <- .Call("call_slra", p, s, r, opt, compute.dp);

  res;    
}


optimize_gsl <- function(sobj, ...) {
  .Call("optimize_gsl", sobj, ...);
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




