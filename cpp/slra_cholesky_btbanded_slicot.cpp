#include <memory.h>
#include <cstdarg>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

#ifdef USE_SLICOT
StationaryCholeskySlicot::
    StationaryCholeskySlicot( const StationaryStructure *s, int D ) :  
      StationaryCholesky(s, D)  {
  myGammaVec = (double*)malloc(getD() * getD() * (getS()+1) * sizeof(double));
  myCholeskyWorkSize = 1 + getS() * getD() * getD() + /* pDW */ 
                       3 * getD() + /* 3 * K */
                       mymax(getS(), getN() - getS()) * getD() * getD();
  myCholeskyWork = (double *)malloc(myCholeskyWorkSize*sizeof(double));                       
}

StationaryCholeskySlicot::~StationaryCholeskySlicot() {
  free(myGammaVec);
  free(myCholeskyWork);
}

void StationaryCholeskySlicot::calcGammaCholesky( const gsl_matrix *R, double reg_gamma  )  {
  size_t info = 0, D = getD(), n = getN();
  const size_t zero = 0;

  computeGammak(R);
  gsl_matrix_vectorize(myGammaVec, myGamma);
    
  mb02gd_("R", "N", &D, &n, &s_minus_1, &zero, 
          &n, myGammaVec, &D, myPackedCholesky, &d_times_s, 
          myCholeskyWork, &myCholeskyWorkSize, &info); /**/

  if (info && reg_gamma > 0) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, 
        "Gamma matrix is singular (MB02GD info = %d), "
        "adding regularization, reg = %f.\n", info, reg_gamma);
    computeGammak(R, reg_gamma);
    gsl_matrix_vectorize(myGammaVec, myGamma);

    mb02gd_("R", "N", &D, &n, &s_minus_1, &zero, 
        &n, myGammaVec, &D, myPackedCholesky, &d_times_s, 
        myCholeskyWork, &myCholeskyWorkSize, &info); /**/
  }

  if (info) {
    throw new Exception("Gamma matrix is singular "
                        "(MB02GD info = %d).\n", info); 
  }
}

#endif  
