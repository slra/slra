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

/*********************/
/* stls: STLS solver */
/*********************/
/* on exit: zero success, otherwise error code: 
   EITER, GSL_ETOLF, GSL_ETOLX, GSL_ETOLG */
/*

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
}*/


static void R_2_x( const gsl_matrix *R, const gsl_matrix *perm, gsl_matrix * x ) {
  gsl_matrix * tempphi = gsl_matrix_alloc(perm->size1, perm->size2);
  gsl_matrix *tempR = gsl_matrix_alloc(R->size1, R->size2);
  size_t minus1 = -1;
  double one = 1;
  double temp;
    
  double *tau = new double[perm->size2];
  size_t status;
  
  gsl_matrix_memcpy(tempphi, perm);
  gsl_matrix_memcpy(tempR, R);
 
 
  PRINTF("R2x entered\n");

  /* Compute LQ factorization of Phi */
  size_t lwork;
  dgelqf_(&(tempphi->size2), &(tempphi->size1), tempphi->data, &(tempphi->tda), 
         tau, &temp, &minus1, &status);
  double *work = new double[(lwork=temp)];
  dgelqf_(&(tempphi->size2), &(tempphi->size1), tempphi->data, &(tempphi->tda), 
         tau, work, &lwork, &status);
  
  delete [] work;


  PRINTF("LQ computed entered\n");
  
  /* Compute R Q^T */
  dormlq_("R", "T", &(tempR->size2), &(tempR->size1), &(tempphi->size2),
         tempphi->data, &(tempphi->size2), tau, tempR->data, &(tempR->tda), 
         &temp, &minus1, &status);
           
  work = new double[(lwork = temp)];
  dormlq_("R", "T", &(tempR->size2), &(tempR->size1), &(tempphi->size2), 
         tempphi->data, &(tempphi->size2), tau, tempR->data, &(tempR->tda), 
         work, &lwork, &status);
  
  gsl_matrix_view RQt = gsl_matrix_submatrix(tempR, 0, 0, tempphi->size2, tempR->size2);
  
  /* solve R_theta L = RQ^T */
  dtrsm_("R", "L", "N", "N", &(RQt.matrix.size2), &(RQt.matrix.size1), &one, 
        tempphi->data,  &(tempphi->tda), RQt.matrix.data, &(RQt.matrix.tda));



  /* Solve AX = B, where  R_theta = [B A] */
  gsl_matrix_view B =  gsl_matrix_submatrix(&(RQt.matrix), 0, 0, 
                           RQt.matrix.size1 - RQt.matrix.size2, RQt.matrix.size2);
  gsl_matrix_view A =  gsl_matrix_submatrix(&(RQt.matrix), RQt.matrix.size1 - RQt.matrix.size2, 0, 
                           RQt.matrix.size2, RQt.matrix.size2);

  size_t *pivot = new size_t[A.matrix.size2];

  /* TODO: Check for singularity of A */
  dgesv_(&(A.matrix.size2), &(B.matrix.size1), A.matrix.data, &(A.matrix.tda), 
         pivot,  B.matrix.data, &(B.matrix.tda), &status);  
         
  gsl_matrix_memcpy(x, &(B.matrix));
  gsl_matrix_scale(x, -1.0);
  
  
  delete [] work;
  delete [] pivot;
  delete [] tau;
  
  gsl_matrix_free(tempphi);
  gsl_matrix_free(tempR);
}


static void compute_lra_R( const gsl_matrix *c, const gsl_matrix *perm, gsl_matrix * R ) {
  size_t status;
  size_t minus1 = -1;
  double temp;

  gsl_matrix * tempc = gsl_matrix_alloc(c->size1, c->size2);
  gsl_matrix_memcpy(tempc, c);
  gsl_matrix * tempu = gsl_matrix_alloc(c->size2, c->size2);
  double *s = new double[mymin(c->size1, c->size2)];
  
  /* Determine optimal work */
  size_t lwork;
  dgesvd_("A", "N", &(tempc->size2), &(tempc->size1),  tempc->data, &(tempc->tda),
         s, tempu->data, &(tempu->size2),NULL, &(tempc->size1), &temp, &minus1, &status);
  double *work = new double[(lwork = temp)];
  /* Compute low-rank approximation */ 
  dgesvd_("A", "N", &(tempc->size2), &(tempc->size1),  tempc->data, &(tempc->tda),
         s, tempu->data, &(tempu->size2), NULL, &(tempc->size1), work, &lwork, &status);

  if (status) {
    delete [] s;  
    delete [] work;  
    gsl_matrix_free(tempc);
    gsl_matrix_free(tempu);

    throw new slraException("Error computing initial approximation: DGESVD didn't converge\n");
  }

  gsl_matrix_transpose(tempu);
  gsl_matrix_view RlraT = gsl_matrix_submatrix(tempu, 0, tempu->size2 - R->size2, tempu->size1, R->size2);
  gsl_matrix_memcpy(R, &(RlraT.matrix));
    
  delete [] s;  
  delete [] work;  
  gsl_matrix_free(tempc);
  gsl_matrix_free(tempu);
}



/* ********************************************************************* * /
/* Interface to the SLICOT function MB02MD. Solves a TLS problem AX = B. * /
/* ********************************************************************* * /
int tls(gsl_matrix* a, gsl_matrix* b, gsl_matrix* x)
{
  double *c, *xa, *s, *dwork, tol=0;
  int    m, n, l, iwarn, info, ldwork, *iwork, ldc, r=0; 

  /* Define constants * /
  m = a->size1;
  n = a->size2;
  l = b->size2;			/* = d * /
  ldc = GSL_MAX(m, n+l);
  if ( m < n + l ) {
    ldwork = GSL_MAX(m * (n + l) + GSL_MAX(3 * m + n + l, 5 * m), 3 * l);
  } else {
    ldwork = GSL_MAX(3*(n+l)+m, 5*(n+l));
  }
  
  /* c = [a b] in FORTRAN format * /
  c = (double *)malloc( ldc * (n+l) * sizeof(double));
  gsl_matrix_vectorize(c, a);
  gsl_matrix_vectorize(c + m*n, b);

  xa = (double *)malloc( n * l * sizeof(double));
  s  = (double *)malloc( (n+l) * sizeof(double));
  dwork = (double *)malloc( ldwork * sizeof(double));
  iwork = (int *)malloc( l * sizeof(int));

  /* Call MB02MD * /
  mb02md_("B", &m, &n, &l, &r, c, &ldc, s, xa, &n, 
	  &tol, iwork, dwork, &ldwork, &iwarn, &info);

  /* copy the result in x * /
  if (!info)
    gsl_matrix_vec_inv(x, xa);

  /* free memory * /
  free(c);
  free(xa);
  free(s);
  free(dwork);
  free(iwork);

  return info;
}


static void checkAndComputeX( slraCostFunction * F, gsl_matrix *x, opt_and_info* opt, gsl_matrix *x_ini ) {
  if (x_ini == NULL) {  /* compute default initial approximation * /
    if (opt->disp == SLRA_OPT_DISP_ITER) {
      PRINTF("X not given, computing TLS initial approximation.\n");
    }
  //  compute_tls2(F->getSMatr(), F->getPerm(), x);
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
}*/

int slra( const gsl_vector *p_in, slraStructure* s, int r, opt_and_info* opt,
         gsl_matrix *x_ini, gsl_matrix *perm, 
         gsl_vector *p_out, gsl_matrix *xh, gsl_matrix *vh ) { 
  slraCostFunction * myCostFun = NULL;
  gsl_matrix *x = NULL;
  gsl_matrix *x2 = NULL;
  gsl_matrix *R = NULL;
  gsl_matrix *U = NULL;
  gsl_rng *rg = NULL;
  int res = GSL_SUCCESS;

  try { 
    myCostFun =  new slraCostFunction(s, r, p_in, opt, perm);
    x = gsl_matrix_alloc(myCostFun->getN(), myCostFun->getD());
    x2 = gsl_matrix_alloc(myCostFun->getN(), myCostFun->getD());
    R = gsl_matrix_alloc(myCostFun->getNplusD(), myCostFun->getD());
    U = gsl_matrix_alloc(myCostFun->getD(), myCostFun->getD());
    
    if (x_ini == NULL) {  /* compute default initial approximation */
      if (opt->disp == SLRA_OPT_DISP_ITER) {
        PRINTF("X not given, computing TLS initial approximation.\n");
      }
      compute_lra_R(myCostFun->getSMatr(), myCostFun->getPerm(), R);
    } else {
      myCostFun->computeR(gsl_matrix_const_submatrix(x_ini, 0, 0, x_ini->size1, x_ini->size2), R);
    }

    R_2_x(R, myCostFun->getPerm(), x);

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
  if (x2 != NULL) {
    gsl_matrix_free(x2);
  }
  if (R != NULL) {
    gsl_matrix_free(R);
  }
  if (rg != NULL) {
    gsl_rng_free(rg);
  }

  if (U != NULL) {
    gsl_matrix_free(U);
  }

  return res;
}

/*#define P ((slra_opt_data*) params)*/


