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

int slra( const gsl_vector *p_in, Structure* s, int r, opt_and_info* opt,
         gsl_matrix *r_ini, gsl_matrix *perm, 
         gsl_vector *p_out, gsl_matrix *rh, gsl_matrix *vh ) { 
  CostFunction * myCostFun = NULL;
  gsl_matrix *x = NULL, *R = NULL;
  int res = GSL_SUCCESS;

  if (opt->gcd) { 
    opt->ls_correction = 1; /* Only correction LS is allowed here */
    if (r_ini == NULL) {
      throw new Exception("GCD computation must have "
                          "an initial approximation.\n");
    }
  }
  
  try { 
    myCostFun =  new CostFunction(s, r, p_in, opt, perm);
    R = gsl_matrix_alloc(myCostFun->getRsize(), myCostFun->getD());
    
    if (r_ini == NULL) {  /* compute default initial approximation */
      myCostFun->computeDefaultRTheta(R);
    } else {
      gsl_matrix_memcpy(R, r_ini);
    }

    x = gsl_matrix_alloc(myCostFun->getRank(), myCostFun->getD());
    myCostFun->R_2_x(R, x);
    time_t t_b = clock();
    gsl_vector_view x_vec = gsl_vector_view_array(x->data, x->size1*x->size2);
    int status = gsl_optimize(myCostFun, opt, &(x_vec.vector), vh);
    opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;
    myCostFun->computeRTheta(x, R);
    if (p_out != NULL) {
      if (opt->gcd) {
        gsl_vector_set_zero(p_out);
      } else {
        if (p_out != p_in) {
          gsl_vector_memcpy(p_out, p_in);
        }
      }
      myCostFun->computeCorrection(p_out, &(x_vec.vector));
    }
    if (rh != NULL) {
      gsl_matrix_memcpy(rh, R);
    }

    gsl_matrix_free(x);
    gsl_matrix_free(R);
    delete myCostFun;
  } catch ( Exception *e ) {
    if (myCostFun != NULL) {
      delete myCostFun;
    }
    if (x != NULL) {
      gsl_matrix_free(x);
    }
    if (R != NULL) {
      gsl_matrix_free(R);
    }
    
    throw; 
  }

  return res;
}
