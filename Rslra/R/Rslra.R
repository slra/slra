 
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
slra <- function(p, s, r = dim(phi)[1] - 1, opt = list(), 
            compute.ph = FALSE, compute.Rh = TRUE, ret.obj = FALSE) {
  # Check necessary parameters            
  if (!is.list(s) || is.null(s$m)) {
    stop('Structure must be a list with "m" and "n" elements');
  } 
  storage.mode(s$m)  <- 'integer';
  storage.mode(s$m)  <- 'double';
  s$m <- as.vector(s$m);
  if (is.null(s$n)) {
    s$n = (length(p) - sum(s$m)) / length(s$m) + 1;
  }
  storage.mode(s$n) <- 'integer';
  storage.mode(s$n) <- 'double';

  if (!prod(s$m > 0) || !prod(s$n > 0)) {
    stop('s$m and s$n elements should be positive vectors');
  }
  
  if (sum(s$m) > sum(s$n)) {
    stop('Matrix should be fat');
  }
  np <- sum(s$m - 1) * length(s$n) + length(s$m) * sum(s$n);
  if (np > length(p)) {
    stop('Vector p too short');  
  } else if (np < length(p)) {
    p <- p[1:np];
    warning('Size of vector p exceeds structure requirements');  
  }
  
  if (is.null(s$phi)) {
    s$phi <- diag(sum(s$m));
  } else{
    if (dim(phi)[1] > dim(phi)[2] || dim(phi)[2] != sum(s$m)) {
      stop('s$phi should be a full row rank matrix compatible with s');
    }
  }
  storage.mode(s$phi) <- 'double';

  if (!is.null(s$w)) {
    if (!prod(s$w > 0)) {
      stop('Weights can be positive or infinite');
    }
  } else {
    s$w <- rep(1, length(s$m));
  }
  storage.mode(s$w) <- 'double';

  storage.mode(r) <- 'integer';
  if (r < 0 || r >= sum(s$m)) {
    stop ('Incorrect r value');
  }
  
  storage.mode(compute.ph) <- storage.mode(compute.Rh) <- 'integer';
  res <- .Call("call_slra", p, s, r, opt, compute.ph, compute.Rh);
}


#optimize_gsl <- function(sobj, ...) {
#  .Call("optimize_gsl", sobj, ...);
#}

#evaluate_slra_f <- function(sobj, ...) {
#  .Call("evaluate_slra_f", sobj, ...);
#}

#evaluate_slra_df <- function(sobj, ...) {
#  .Call("evaluate_slra_df", sobj, ...);
#}

#return_slra_object <- function(sobj) {
#  .Call("return_slra_object", sobj);
#}

#compute_slra_correction <- function(sobj) {
#  .Call("compute_slra_correction", sobj);
#}/*


