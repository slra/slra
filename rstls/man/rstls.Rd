\name{rstls}
\alias{stls}

\title{Function for solving structured total least squares problem}

\description{R interface for the STLS package solver}

\usage{
stls <- function(A = NULL, B = NULL, S, X = NULL, opts = list(epsabs = 0, epsrel = 1e-6, epsgrad = 1e-6, maxiter = 100, disp = 'iter'),
           P = NULL, compute.dp = !is.null(P))
}

\arguments{
  \item{A}{A matrix}
  \item{B}{B matrix}
  \item{S}{structure specification}
  \item{X}{initial approximation}
  \item{opts}{Options for underlying local optimization algorithm}
  \item{P}{vector of parameters (e.g. time series for Hankel matrices) }
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
    \item{m}{-- total number of rows. Can be omitted if the structured matrix is given by (A,B)}
    \item{k}{-- the size of sub blocks in each 'T'/'H' block}
    \item{A}{-- A q\times 3 array (q <= 10), where each row describes a block
      A[i,1] = 1 for 'T', 2 for 'H', 3 for 'U', 4 for 'E'
      A[i,2] = Total number of columns in a block
      A[i,3] = Number of columns in subblocks (if block is 'T' or 'H')
    }
    \item{d}{The rank reduction. Can be omitted if either (A,B) or X is given
  }
  
  Note: the STLS approximation is computed only if the structured matrix was given by vector P.
}


\example{
# Computing STLS approximation of a Hankel matrix with default initial approximation

f <- sin(1:100 * (2 * pi /10)) + norm(length(f0), 0, 0.1;
res <- stls(P = f, S=list(m=(length(f)-r), k=1, d=1, A=c(2, r+1, 1)), X = xini);
print(res);

}



\references{
  I. Markovsky and S. {Van Huffel}
  Software for structured total least squares estimation: User's guide
  Technical report 03--136, Dept. EE, K.U. Leuven, 2002.
}
