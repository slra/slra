\name{Rslra}
\alias{slra}

\title{Function for solving structured low rank approximation problem}

\description{R interface for the SLRA package solver}

\usage{
stls <- function(P, S, R, X = NULL, opts = list(epsabs = 0, epsrel = 1e-5, epsgrad = 1e-5, maxiter = 100, disp = 'iter'),
           compute.dp = FALSE)
}

\arguments{
  \item{P}{vector of parameters (e.g. time series for Hankel matrices) }
  \item{S}{structure specification}
  \item{R}{rank (by default the number of columns - 1)}
  \item{X}{initial approximation}
  \item{opts}{Options for underlying local optimization algorithm}
  \item{compute.dp}{whether to compute the approximation obtained by STLS}
}

\value{
  The returned value is a list with components
  \item{xh}{A matrix containing the local minimum point}
  \item{info}{Output information  as describe in (Markovsky, Van Huffel, 2002)}
  \item{vh}{The covariance matrix for the X estimate}
  \item{ph}{The STLS estimate of the vector P.}
}

\details{
  The core C function and the package for STLS problem is described in (Markovsky, Van Huffel, 2002), please read it first.
  
  This interface has additional features (compared to (Markovsky, Van Huffel, 2002)):
   however 
  \describe{
    \item{1}{the structured matrix can be given not only by A and B but also by p}
    \item{2}{the stls approximation can be computed (not only optimal X as before)}
  }
  
  
  The structure specification should be currently a list with following elements:
  \describe{
    \item{k}{-- the size of sub blocks in each 'T'/'H' block}
    \item{A}{-- A q\times 3 array (q <= 10), where each row describes a block
      A[i,1] = 1 for 'T', 2 for 'H', 3 for 'U', 4 for 'E'
      A[i,2] = Total number of columns in a block
      A[i,3] = Number of columns in subblocks (if block is 'T' or 'H')
    }
  }
  If k = 1, then the structure can be just a matrix (or a vector).
  
  Note: the STLS approximation is computed only if the structured matrix was given by vector P.
}


\examples{
# Computing STLS approximation of a Hankel matrix with default initial approximation
library('Rslra');
r <- 2;
f <- sin(1:100 * (2 * pi /10)) + 0.03 * rnorm(100, 0, 0.1);
res <- slra(f, c(2, r+1, 1), r, compute.dp = TRUE);
print(res);

}



\references{
  I. Markovsky and S. {Van Huffel}
  Software for structured total least squares estimation: User's guide
  Technical report 03--136, Dept. EE, K.U. Leuven, 2002.
}
