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
  CostFunction * myCostFun = NULL;
  gsl_matrix *x = NULL, *Rtheta = NULL, *tmpphi = NULL;
  
  if (Psi != NULL) {
    if (Phi != NULL) {
      tmpphi = gsl_matrix_alloc(Phi->size1, Psi->size2);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Phi, Psi, 0, tmpphi); 
    } else {
      tmpphi = gsl_matrix_alloc(Psi->size1, Psi->size2);
      gsl_matrix_memcpy(tmpphi, Psi); 
    }
  }

  if (opt->gcd) { 
    opt->ls_correction = 1; /* Only correction LS is allowed for GCD */
    if (Rini == NULL) {
      throw new Exception("GCD computation must have "
                          "an initial approximation.\n");
    }
  }
  
  try { 
    time_t t_b = clock();
    myCostFun =  new CostFunction(s, d, p_in, opt, 
                                  (tmpphi != NULL ? tmpphi : Phi));
    Rtheta = gsl_matrix_alloc(myCostFun->getRsize(), myCostFun->getD());

    if (Rini == NULL) {  /* compute default initial approximation */
      myCostFun->computeDefaultRTheta(Rtheta);
    } else {
      if (tmpphi == NULL) {
        gsl_matrix_memcpy(Rtheta, Rini);
      } else {
        ls_solve(Psi, Rini, Rtheta);
      }      
    }

    x = gsl_matrix_alloc(myCostFun->getRank(), myCostFun->getD());
    myCostFun->Rtheta2X(Rtheta, x);
    gsl_vector_view x_vec = gsl_vector_view_array(x->data, x->size1*x->size2);
    int status = gsl_optimize(myCostFun, opt, &(x_vec.vector), vh);
    if (p_out != NULL) {
      if (p_out != p_in && !opt->gcd) {
        gsl_vector_memcpy(p_out, p_in);
      }
      myCostFun->computeCorrection(p_out, &(x_vec.vector));
    }

    opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;

    if (Rout != NULL) {
      myCostFun->X2Rtheta(x, Rtheta);
      if (Psi == NULL) {
        gsl_matrix_memcpy(Rout, Rtheta);
      } else {
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Psi, Rtheta, 0, Rout);
      }
    }

    throw (Exception *)NULL; /* Throw NULL exception to unify deallocation */
  } catch ( Exception *e ) {
    if (myCostFun != NULL) {
      delete myCostFun;
    }
    if (x != NULL) {
      gsl_matrix_free(x);
    }
    if (Rtheta != NULL) {
      gsl_matrix_free(Rtheta);
    }
    if (tmpphi != NULL) {
      gsl_matrix_free(tmpphi);
    }
    
    if (e != NULL) { /* Abnormal termination only if e is normal exception */
      throw;  
    }
  }
}
