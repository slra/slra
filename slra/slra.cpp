/* stls.c: implementations of the functions from stls.h */

#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>


#include "slra.h"

/*********************/
/* stls: STLS solver */
/*********************/
/* on exit: zero success, otherwise error code: 
   EITER, GSL_ETOLF, GSL_ETOLX, GSL_ETOLG */


static int compute_tls( const gsl_matrix *c, const gsl_matrix *perm, gsl_matrix * x ) {
  gsl_matrix * tempc = gsl_matrix_alloc(c->size1, c->size2);
  gsl_matrix_view c_sub_a, c_sub_b;
  int n = x->size1, d = x->size2;
  int status;
    
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, c, perm, 0.0, tempc);
  c_sub_a = gsl_matrix_submatrix(tempc, 0, 0, c->size1, n);
  c_sub_b = gsl_matrix_submatrix(tempc, 0, n, c->size1, d);
    
  status = tls(&c_sub_a.matrix, &c_sub_b.matrix, x);
    
  gsl_matrix_free(tempc);

  if (status) {
    PRINTF("Initial approximation can not be computed: "
        "MB02MD failed with an error info = %d.\n", status);

    return GSL_EINVAL;
  }

  return GSL_SUCCESS;
}


/* ********************************************************************* */
/* Interface to the SLICOT function MB02MD. Solves a TLS problem AX = B. */
/* ********************************************************************* */
int tls(gsl_matrix* a, gsl_matrix* b, gsl_matrix* x)
{
  double *c, *xa, *s, *dwork, tol=0;
  int    m, n, l, iwarn, info, ldwork, *iwork, ldc, r=0; 

  /* Define constants */
  m = a->size1;
  n = a->size2;
  l = b->size2;			/* = d */
  ldc = GSL_MAX(m, n+l);
  if ( m < n + l ) {
    ldwork = GSL_MAX(m * (n + l) + GSL_MAX(3 * m + n + l, 5 * m), 3 * l);
  } else {
    ldwork = GSL_MAX(3*(n+l)+m, 5*(n+l));
  }
  
  /* c = [a b] in FORTRAN format */
  c = (double *)malloc( ldc * (n+l) * sizeof(double));
  gsl_matrix_vectorize(c, a);
  gsl_matrix_vectorize(c + m*n, b);

  xa = (double *)malloc( n * l * sizeof(double));
  s  = (double *)malloc( (n+l) * sizeof(double));
  dwork = (double *)malloc( ldwork * sizeof(double));
  iwork = (int *)malloc( l * sizeof(int));

  /* Call MB02MD */
  mb02md_("B", &m, &n, &l, &r, c, &ldc, s, xa, &n, 
	  &tol, iwork, dwork, &ldwork, &iwarn, &info);

  /* copy the result in x */
  if (!info)
    gsl_matrix_vec_inv(x, xa);

  /* free memory */
  free(c);
  free(xa);
  free(s);
  free(dwork);
  free(iwork);

  return info;
}


static void checkAndComputeX( slraCostFunction * F, gsl_matrix *x, opt_and_info* opt, gsl_matrix *x_ini ) {
  if (x_ini == NULL) {  /* compute default initial approximation */
    if (opt->disp == SLRA_OPT_DISP_ITER) {
      PRINTF("X not given, computing TLS initial approximation.\n");
    }
    if (compute_tls(F->getSMatr(), F->getPerm(), x) ==  GSL_EINVAL) {
      throw new slraException("Error while computing initial approximation.\n");   
    }
  } else {
    if (x_ini->size1 != x->size1 || x_ini->size2 != x->size2) {
      throw new slraException("Initial approximation doesn't conform to " 
                            "the structure specification.\n");   
    }      
    gsl_matrix_memcpy(x, x_ini);
  }
}

int slra( const gsl_vector *p_in, slraStructure* s, int r, opt_and_info* opt,
         gsl_matrix *x_ini, gsl_matrix *perm, 
         gsl_vector *p_out, gsl_matrix *xh, gsl_matrix *vh ) { 
  slraCostFunction * myCostFun = NULL;
  gsl_matrix *x = NULL;
  int res = GSL_SUCCESS;

  try { 
    myCostFun =  new slraCostFunction(s, r, p_in, opt, perm);
    x = gsl_matrix_alloc(myCostFun->getN(), myCostFun->getD());
    
    checkAndComputeX(myCostFun, x, opt, x_ini);
    
    time_t t_b = clock();
    gsl_vector_view x_vec = gsl_vector_view_array(x->data, x->size1 * x->size2);
    int status = slra_gsl_optimize(myCostFun, opt, &(x_vec.vector), vh);
    opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;

    if (p_out != NULL) {
      if (p_out != p_in) {
        gsl_vector_memcpy(p_out, p_in);
      }
      myCostFun->computeCorrection(p_out, &(x_vec.vector));
    }
    if (xh != NULL) {
      gsl_matrix_memcpy(xh, x);
    }
  } catch ( slraException *e ) {
    res = GSL_EINVAL;
    PRINTF(e->getMessage());
    delete e;
  }

  if (myCostFun != NULL) {
    delete myCostFun;
  }
  
  if (x != NULL) {
    gsl_matrix_free(x);
  }

  return res;
}

/*#define P ((slra_opt_data*) params)*/


