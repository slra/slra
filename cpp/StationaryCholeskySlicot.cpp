#include <memory.h>
#include <cstdarg>
#include "slra.h"

#ifdef USE_SLICOT
StationaryCholeskySlicot::
    StationaryCholeskySlicot( const StationaryStructure *s, size_t d ) :  
      StationaryCholesky(s, d)  {
  myGammaVec = (double*)malloc(d * d * (getMu() + 1) * sizeof(double));
  myCholeskyWorkSize = 1 + getMu() * d * d + /* pDW */ 
                       3 * d + /* 3 * K */
                       mymax(getMu(), getN() - getMu()) * d * d;
  myCholeskyWork = (double *)malloc(myCholeskyWorkSize * sizeof(double));
}

StationaryCholeskySlicot::~StationaryCholeskySlicot() {
  free(myGammaVec);
  free(myCholeskyWork);
}

void StationaryCholeskySlicot::calcGammaCholesky( const gsl_matrix *Rt, double reg  )  {
  size_t info = 0, d = getD(), n = getN();
  const size_t zero = 0;

  computeGammak(Rt);
  gsl_matrix_vectorize(myGammaVec, myGamma);
    
  mb02gd_("R", "N", &d, &n, &myMu_1, &zero, 
          &n, myGammaVec, &d, myPackedCholesky, &myDMu, 
          myCholeskyWork, &myCholeskyWorkSize, &info); /**/

  if (info && reg_gamma > 0) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, 
        "Gamma matrix is singular (MB02GD info = %d), "
        "adding regularization, reg = %f.\n", info, reg_gamma);
    computeGammak(R, reg_gamma);
    gsl_matrix_vectorize(myGammaVec, myGamma);

    mb02gd_("R", "N", &D, &n, &&myMu_1, &zero, 
        &n, myGammaVec, &D, myPackedCholesky, &myDMu, 
        myCholeskyWork, &myCholeskyWorkSize, &info); /**/
  }

  if (info) {
    throw new Exception("Gamma matrix is singular "
                        "(MB02GD info = %d).\n", info); 
  }
}

#endif  
