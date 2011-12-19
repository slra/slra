\name{Rslra}
\alias{slra}

\title{Function for solving structured low rank approximation problem}

\description{SLRA - solves the structured low rank approximation problem 
minimize over Ph  ||P - Ph||_2 subject to rank(S(Ph)) <= r
where the matrix S(Ph) is of the form [S1 ... Sq], with Si being 
(block) Toeplitz, Hankel structured, unstructured, or constant}

\usage{
stls <- function(P, S, R, X0, opt, compute.dp)
}

\arguments{
  \item{P}{parameter vector}
  \item{S}{structure specification}
  \item{R}{rank (default is rank reduction by 1)}
  \item{X0}{initial approximation}
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
  
  The structure specification should be currently a q x 1-4 matrix,
  \describe{
    \item{S(i, 1)}{= number of subblocks in a block row of Si 
                     (use S(i, 1) = 1 for unstructured block Si)}
    \item{S(i, 2)}{= number of columns in a subblock (default 1)}
    \item{S(i, 3)}{= 1 if the block is exact, or 0 otherwise (default 0)}
    \item{S(i, 4)}{= 1 if the block is Toeplitz, 0 if it is Hankel (default 0)}
  }
  For block-Toeplitz/Hankel structured block Si with K x nu(i) size blocks
  S is a list with fields k and a, where S$A is the matrix described 
  above and S$k = K (the same for all block-Toeplitz/Hankel blocks S1)

  opt is a list with elements:
  \describe{
    \item{Most widely used:}{}
    \item{OPT.DISP}{- information about progress of the optimization 
               options: 'iter', 'notify', 'final', 'off' (default 'notify')}
    \item{OPT.MAXITER}{- maximum number of iterations,}
    \item{OPT.METHOD}{- optimization method (consult GSL library manual): 
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
    \item{Initial parameters of optimization methods:}{}
    \item{OPT.STEP}{'step_size' for fdfminimizer_set, fminimizer_set (scalar)} 
    \item{OPT.TOL}{'tol' for fdfminimizer_set, fminimizer_set,} 
    \item{Stopping parameters:}
    \item{OPT.EPSREL, OPT.EPSABS}{- 'gsl_multifit_test_delta' stopping criterion}
    \item{OPT.EPSGRAD}{- 'gsl_multi..._test_gradient' stopping criterion}
    \item{Advanced parameters:}{}
    \item{OPT.REGGAMMA}{ - regularization parameter for gamma, absolute}
  }      
  
  Note: the STLS approximation is computed only if the structured matrix was given by vector P.
}


\examples{
# Computing STLS approximation of a Hankel matrix with default initial approximation
library('Rslra');
r <- 2;
f <- sin(1:100 * (2 * pi /10)) + 0.03 * rnorm(100, 0, 0.1);
res <- slra(f, r+1, compute.dp = TRUE);
print(res);

}



\references{
  I. Markovsky and S. {Van Huffel}
  Software for structured total least squares estimation: User's guide
  Technical report 03--136, Dept. EE, K.U. Leuven, 2002.
}
