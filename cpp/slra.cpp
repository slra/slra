/* stls.c: implementations of the functions from stls.h */

#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "slra.h"

void slra( const gsl_vector *p_in, Structure* s, int d, 
          OptimizationOptions* opt, gsl_matrix *Rini, gsl_matrix *Phi, 
          gsl_matrix *Psi, gsl_vector *p_out, gsl_matrix *Rout, 
          gsl_matrix *vh ) { 
  OptFunctionSLRA *optFun = NULL;
  CostFunction * myCostFun = NULL;
  gsl_vector *x = NULL;
  
  try { 
    time_t t_b = clock();
    myCostFun =  new CostFunction(s, d, p_in, opt->reggamma);

    if ( opt->ls_correction) {
      optFun = new OptFunctionSLRACorrection(*myCostFun, Phi, Psi);
    } else { 
      optFun = new OptFunctionSLRACholesky(*myCostFun, Phi, Psi);
    }

    x = gsl_vector_alloc(optFun->getNvar());

    if (Rini == NULL) {  
      Log::lprintf(Log::LOG_LEVEL_ITER, 
           "R not given - computing initial approximation.\n");    
      optFun->computeDefaultx(x);
    } else {
      optFun->RTheta2x(Rini, x);
    }

    gsl_optimize(optFun, opt, x, vh);
    if (p_out != NULL) {
      if (p_out != p_in) {
        gsl_vector_memcpy(p_out, p_in);
      }
      optFun->computeCorrection(p_out, x);
    }

    opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;
    
    if (Rout != NULL) {
      optFun->x2RTheta(Rout, x);
    }

    throw (Exception *)NULL; /* Throw NULL exception to unify deallocation */
  } catch ( Exception *e ) {
    if (optFun != NULL)  {
      delete optFun;
    }
    if (myCostFun != NULL) {
      delete myCostFun;
    }
    if (x != NULL) {
      gsl_vector_free(x);
    }
    
    if (e != NULL) { /* Abnormal termination only if e is normal exception */
      throw;  
    }
  }
}
