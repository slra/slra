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


void slra( VarproFunction *costFun, 
           OptimizationOptions* opt, gsl_matrix *Rini, gsl_matrix *Psi, 
           gsl_vector *p_out, gsl_matrix *r_out, gsl_matrix *v_out ) { 
  OptFunctionSLRA *optFun = NULL;
  gsl_vector *x = NULL;
  double old_reg = costFun->getReggamma();
  
  try { 
    time_t t_b = clock();

    costFun->setReggamma(opt->reggamma);
    if (opt->ls_correction || costFun->isGCD()) {
      optFun = new OptFunctionSLRACorrection(*costFun, Psi);
    } else { 
      optFun = new OptFunctionSLRACholesky(*costFun, Psi);
    }

    x = gsl_vector_alloc(optFun->getNvar());

    if (Rini == NULL) {  
      Log::lprintf(Log::LOG_LEVEL_ITER, 
           "R not given - computing initial approximation.\n");    
      optFun->computeDefaultx(x);
    } else {
      optFun->RTheta2x(Rini, x);
    }

    opt->gslOptimize(optFun, x, v_out);

    if (p_out != NULL) {
      optFun->computePhat(p_out, x);
    }
    opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;
    if (r_out != NULL) {
      optFun->x2RTheta(r_out, x);
    }

    throw (Exception *)NULL; /* Throw NULL exception to unify deallocation */
  } catch ( Exception *e ) {
    if (optFun != NULL)  {
      delete optFun;
    }
    gsl_vector_free_ifnull(x);
    
    if (e != NULL) { /* Abnormal termination only if e is normal exception */
      throw;  
    }
    costFun->setReggamma(old_reg);
  }
}
