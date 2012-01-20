library('Rslra');
library("minpack.lm");
library('stats');
library('gsl');


slra_f_nls <- function(par, sobj) { evaluate_slra_f(sobj, par); };
slra_df_nls <- function(par, sobj) { evaluate_slra_df(sobj, par); };

slra_f <- function(par, sobj) { 
  fv <- evaluate_slra_f(sobj, par); 
  sum(fv * fv); 
};

slra_df <- function(par, sobj) { 
  fv <- evaluate_slra_f(sobj, par); 
  jv <- evaluate_slra_df(sobj, par); 
  2 * t(jv) %*% fv; 
};

norm2 <- function(f) {
  sqrt(sum(f * f));
}

ext_optimize <- function(sobj, fun, ..., is.nls = FALSE) {
  x0 <- as.vector(return_slra_object(sobj)$xh);
  if (is.nls) {
    fun(par = x0, fn = slra_f_nls, jac = slra_df_nls, sobj = sobj, ...);
  } else {
    fun(par = x0, fn = slra_f, gr = slra_df, sobj = sobj, ...);
  }
}

optimize_gsl <- function(par, sobj, ...) {
  state <-  multimin.init(x=par, f = function(x) { slra_f(x, sobj); }, df=function(x) { slra_df(x, sobj); }, ...); 
  
  for (i in 1:100) {
    state <- multimin.iterate(state);
    grnorm <- norm2(state$df);
    print(c(i-1, state$f, grnorm)); 
    if (grnorm < 1e-4) {
      break;
    }
  }
  state;
}


f0 <- 1:40;
r <- 2;
f <- f0 + rnorm(length(f0)) * 3;

slralmres <- slra(f, r+1, opt=list(disp="iter"));
sobj <- slra(f, r+1, opt=list(disp="iter"), ret.obj=TRUE);
x0 <- as.vector(return_slra_object(sobj)$xh);
lmres <- nls.lm(par = x0, fn = slra_f_nls, jac = slra_df_nls, sobj = sobj, control=nls.lm.control(nprint=1));


slraqnres <- slra(f, r+1, opt=list(disp="iter",method="q2"));

sobj <- slra(f, r+1, opt=list(disp="iter"), ret.obj=TRUE);

x0 <- as.vector(return_slra_object(sobj)$xh);
optimres <- optim(par = x0, fn = slra_f, gr = slra_df, sobj = sobj, method="BFGS", control=list(REPORT=1,trace=3));
fls2 <- evaluate_slra_f(sobj, optimres$xh);
compute_slra_correction(sobj)


sobj <- slra(f, r+1, opt=list(disp="iter"), ret.obj=TRUE);
gslres <- optimize_gsl(par = x0, sobj = sobj, method="bfgs", tol=1e-6, step.size=0.001);



#sobj <- slra(f, r+1, opt=list(disp="iter"), ret.obj=TRUE);
#optimize_slra(sobj);
#res <- return_slra_object(sobj);
#fls2 <- evaluate_slra_f(sobj, res$xh);
#print(sum(fls2*fls2));

#print(return_slra_object(sobj));
#print(compute_slra_correction(sobj));


