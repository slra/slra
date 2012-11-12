#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

OptFunctionSLRA::OptFunctionSLRA( CostFunction &fun, gsl_matrix *perm ) : myFun(fun) {
  if (perm == NULL) {
    myPerm = gsl_matrix_alloc(myFun.getM(), myFun.getM());
    gsl_matrix_set_identity(myPerm);
  } else {
    if (perm->size1 != myFun.getM() || perm->size1 < perm->size2) {
      throw new Exception("Incorrect sizes of permutation matrix.\n");   
    }
    myPerm = gsl_matrix_alloc(perm->size1, perm->size2);
    gsl_matrix_memcpy(myPerm, perm);
  }
  myTmpR = gsl_matrix_alloc(myFun.getM(), myFun.getD());
  myTmpGradR = gsl_matrix_alloc(myFun.getM(), myFun.getD());
  myTmpThetaExt = gsl_matrix_alloc(myPerm->size2, myFun.getD());
}

OptFunctionSLRA::~OptFunctionSLRA()  {
  gsl_matrix_free(myTmpR);
  gsl_matrix_free(myTmpGradR);
  gsl_matrix_free(myPerm);
  gsl_matrix_free(myTmpThetaExt);
}

int OptFunctionSLRA::getNvar() { 
  return myFun.getRank() * myFun.getD(); 
}

void OptFunctionSLRA::X2Rtheta( const gsl_matrix *x, gsl_matrix *RTheta ) { 
  gsl_vector diag;
  gsl_matrix submat;
  int n = x->size1, d = x->size2;

  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  submat = gsl_matrix_submatrix(RTheta, 0, 0, n, d).matrix;
  gsl_matrix_memcpy(&submat, x); 
  submat = gsl_matrix_submatrix(RTheta, n, 0, d, d).matrix;
  gsl_matrix_set_all(&submat, 0);
  diag   = gsl_matrix_diagonal(&submat).vector;    
  gsl_vector_set_all(&diag, -1);
}

void OptFunctionSLRA::computeR( const gsl_vector * x, gsl_matrix *R ) { 
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector(x, myFun.getRank(), myFun.getD());

  X2Rtheta(&x_mat.matrix, myTmpThetaExt);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPerm, myTmpThetaExt, 0, R);
}

void OptFunctionSLRA::computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad ) {
  computeR(x, myTmpR);

  if (grad == NULL) {
    myFun.computeFuncAndGrad(myTmpR, f, NULL);
  } else {
    myFun.computeFuncAndGrad(myTmpR, f, myTmpGradR);

    gsl_matrix_view grad_matr = gsl_matrix_view_vector(grad, myFun.getRank(), myFun.getD());
    gsl_matrix_view perm_sub_matr = gsl_matrix_submatrix(myPerm, 0, 0, 
                                      myFun.getM(), myFun.getRank());
    /* Compute gradient of f_{\Phi}(X) */ 
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &perm_sub_matr.matrix, 
                   myTmpGradR, 0.0, &grad_matr.matrix);
  }                   
}   

