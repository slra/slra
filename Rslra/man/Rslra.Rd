\name{Rslra}
\alias{slra}

\title{Function for solving structured low rank approximation problem}

\description{SLRA - solves the structured low rank approximation problem 

\deqn{\mbox{minimize over}\;\widehat{p}\quad\|p - \widehat{p}\|_w\quad
\mbox{subject to}\;\mathrm{rank}(\mathcal{S}(\widehat{p})) \le r}{minimize 
over ph ||p - ph||_w subject to  rank S(p) <= r}

where \eqn{\mathcal{S}(p)}{S(p)} is the affine structure of type 
\eqn{\Phi \mathcal{H}_{m,n}}{\Phi H_mn}, 
where \eqn{H_{m,n}} is a Mosaic Hankel matrix.
}

\usage{
stls <- function(p, s, r, opt, compute.ph, compute.Rh)
}

\arguments{
  \item{p}{parameter vector of length \eqn{n_p}}
  \item{s}{structure specification \eqn{\mathcal{S}(p)}}{S(p)}
  \item{r}{rank (default is rank reduction by 1)}
  \item{opt}{optimization parameters}
  \item{compute.dp}{whether to compute the approximation obtained by STLS}
}

\value{
  The returned value is a list with components
  \item{xh}{low rank certificate: S(Ph) * [XH; -I] = 0}
  \item{info}{convergence information: 
    info$iter - number of iterations
    info$time - execution time
    info$fmin - (||P - Ph||_2)^2}
  \item{vh}{Asymptotic covariance matrix of XH}
  \item{ph}{Ph vector (approximation)}
}

\details{
  The core C function and the package for STLS problem is described in 
  (Markovsky, Van Huffel, 2002), please read it first.
  
  The structure specification \code{s} should be a list with following fields
  \describe{
    \item{m}{(required) vector m of length q for the Mosaic structure}
    \item{n}{- vector n of length N for the Mosaic structure}
    \item{phi}{- (optional) matrix \eqn{\Phi} (default \eqn{I_d})}
    \item{w}{- (optional) weights. Can be:
       \describe{
        \item{vector of length n_p}{elementwise weights}
        \item{q x N matrix}{weights for each block}
        \item{vector of length q}{weights for each block row}
      }
    }
  }

  Optimization parameters \code{opt} are passed in a list:
  \describe{
    \item{Most widely used:}{
      \describe{
    \item{disp}{- information about progress of the optimization 
               options: 'iter', 'notify', 'final', 'off' (default 'notify')}
    \item{maxiter}{- maximum number of iterations,}
    \item{method}{- optimization method (consult GSL library manual): 
      \describe{
        \item{'l'}{Levenberg-Marquardt, gsl_multifit_fdf_solver_...}
        \item{'ll'}{...lmder (default)}
        \item{'ls'}{...lmsder}
        \item{'q'}{Quasi-Newton, gsl_multimin_fdfminimizer_...}
        \item{'qb'}{...bfgs (default)}
        \item{'q2'}{...bfgs2}
        \item{'qp'}{...conjugate_pr}
        \item{'qf'}{...conjugate_fr}
        \item{'n'}{Nead-Melder, gsl_multimin_fminimizer_...}
        \item{'nn'}{...nmsimplex (default)}
        \item{'n2'}{...nmsimplex2}
        \item{'nr'}{...nmsimplex2rand}
      }
    }        
    }}
    \item{Initial parameters of optimization methods:}{}
    \item{step}{'step_size' for fdfminimizer_set, fminimizer_set (scalar)} 
    \item{tol}{'tol' for fdfminimizer_set, fminimizer_set,} 
    \item{Stopping parameters:}{}
    \item{epsrel,opt$epsabs}{- 'gsl_multifit_test_delta' stopping criterion}
    \item{epsgrad}{- 'gsl_multi..._test_gradient' stopping criterion}
    \item{Advanced parameters:}{}
    \item{reggamma}{ - regularization parameter for gamma, absolute}
  }      
  
  Note: the STLS approximation is computed only if the structured matrix was given by vector P.
}


\examples{
# Computing STLS approximation of a Hankel matrix with default initial approximation
library('Rslra');
r <- 2;
T <- 100;
f <- sin(1:T * (2 * pi /10)) + 0.03 * rnorm(T, 0, 0.1);
res <- slra(f, list(m = r+1, n = T-r), r, compute.ph = TRUE);
print(res);
}



\references{
  I. Markovsky and S. {Van Huffel}
  Software for structured total least squares estimation: User's guide
  Technical report 03--136, Dept. EE, K.U. Leuven, 2002.
}
