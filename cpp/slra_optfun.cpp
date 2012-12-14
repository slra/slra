#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

OptFunctionSLRA::OptFunctionSLRA( CostFunction &fun, gsl_matrix *Psi ) : 
    myFun(fun) {
  if (Psi == NULL) {
    myPsi = gsl_matrix_alloc(myFun.getNrow(), myFun.getNrow());
    gsl_matrix_set_identity(myPsi);
  } else {
    myPsi = gsl_matrix_alloc(Psi->size1, Psi->size2);
    if (myPsi->size1 != myFun.getNrow() || myPsi->size1 < myPsi->size2) {
      throw new Exception("Incorrect sizes of Psi matrix.\n");   
    }
    gsl_matrix_memcpy(myPsi, Psi); 
  }
  
  myTmpR = gsl_matrix_alloc(myFun.getNrow(), myFun.getD());
  myTmpXId = gsl_matrix_alloc(myPsi->size2, myFun.getD());
}

OptFunctionSLRA::~OptFunctionSLRA()  {
  gsl_matrix_free(myTmpR);
  gsl_matrix_free(myTmpXId);
  gsl_matrix_free(myPsi);
}

size_t OptFunctionSLRA::getNvar() { 
  return getRank() * myFun.getD(); 
}

gsl_matrix OptFunctionSLRA::x2xmat( const gsl_vector *x ) {
  return gsl_matrix_const_view_vector(x, getRank(), myFun.getD()).matrix;
}

void OptFunctionSLRA::computeR( const gsl_vector * x, gsl_matrix *R ) { 
  gsl_matrix x_mat = x2xmat(x);
  X2XId(&x_mat, myTmpXId);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPsi, myTmpXId, 0, R);
}

void OptFunctionSLRA::computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad ) {
  computeR(x, myTmpR);
  if (grad == NULL) {
    myFun.computeFuncAndGrad(myTmpR, f, NULL, NULL);
  } else {
    gsl_matrix grad_matr = gsl_matrix_view_vector(grad, getRank(), myFun.getD()).matrix;
    gsl_matrix psi_sub_matr = gsl_matrix_submatrix(myPsi, 0, 0, 
                                      myFun.getNrow(), getRank()).matrix;
    myFun.computeFuncAndGrad(myTmpR, f, &psi_sub_matr, &grad_matr);
  }                   
}   

void OptFunctionSLRACholesky::computeFuncAndJac( const gsl_vector* x, 
         gsl_vector *res, gsl_matrix *jac ) {
  gsl_matrix psi_sub_matr = gsl_matrix_submatrix(myPsi, 0, 0, 
                                 myFun.getNrow(), getRank()).matrix;         
  computeR(x, myTmpR);
  myFun.computeFuncAndPseudoJacobianLs(myTmpR, &psi_sub_matr, res, jac); 
}   

void OptFunctionSLRACorrection::computeFuncAndJac( const gsl_vector* x, 
         gsl_vector *res, gsl_matrix *jac ) {
  gsl_matrix psi_sub_matr = gsl_matrix_submatrix(myPsi, 0, 0, 
                                 myFun.getNrow(), getRank()).matrix;         
  computeR(x, myTmpR);
  myFun.computeCorrectionAndJacobian(myTmpR, &psi_sub_matr, res, jac); 
}   

void OptFunctionSLRA::X2XId( const gsl_matrix *x, gsl_matrix *XId ) { 
  gsl_vector diag;
  gsl_matrix sm;
  size_t n = x->size1, d = x->size2;

  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  gsl_matrix_memcpy(&(sm = gsl_matrix_submatrix(XId, 0, 0, n, d).matrix), x); 
  gsl_matrix_set_all(&(sm = gsl_matrix_submatrix(XId, n, 0, d, d).matrix), 0);
  gsl_vector_set_all(&(diag = gsl_matrix_diagonal(&sm).vector), -1);
}

void OptFunctionSLRA::PQ2XId( const gsl_matrix *R, gsl_matrix * x ) {
  gsl_matrix *tR = gsl_matrix_alloc(R->size1, R->size2);
  gsl_matrix_memcpy(tR, R);
  size_t status = 0, s1_s2 = tR->size1 - tR->size2;
  /* Solve AX = B, where  R_theta = [B A] */
  gsl_matrix B = gsl_matrix_submatrix(tR, 0, 0, s1_s2, tR->size2).matrix;
  gsl_matrix A = gsl_matrix_submatrix(tR, s1_s2, 0, tR->size2, tR->size2).matrix;
  size_t *pivot = new size_t[A.size2];
  /* TODO: Check for singularity of A */
  dgesv_(&A.size2, &B.size1, A.data, &A.tda, pivot, B.data, &B.tda, &status);  
  delete [] pivot;
  gsl_matrix_memcpy(x, &B);
  gsl_matrix_scale(x, -1.0);
  gsl_matrix_free(tR);
}

void OptFunctionSLRA::RTheta2x( gsl_matrix *RTheta, gsl_vector *x ) {
  gsl_matrix x_mat = x2xmat(x);
  if (myPsi == NULL) {
    gsl_matrix_memcpy(myTmpXId, RTheta);
  } else {
    ls_solve(myPsi, RTheta, myTmpXId);
  }
  PQ2XId(myTmpXId, &x_mat);
}

void OptFunctionSLRA::x2RTheta( gsl_matrix *RTheta, const gsl_vector *x ) {
  gsl_matrix x_mat = x2xmat(x);
  X2XId(&x_mat, myTmpXId);
  if (myPsi == NULL) {
    gsl_matrix_memcpy(RTheta, myTmpXId);
  } else {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPsi, myTmpXId, 0, RTheta);
  }
}

void OptFunctionSLRA::computeDefaultx( gsl_vector *x ) {
 gsl_matrix *Rtheta = gsl_matrix_alloc(myFun.getNrow(), myFun.getD());
 myFun.computeDefaultRTheta(Rtheta);
 RTheta2x(Rtheta, x);
 gsl_matrix_free(Rtheta);
}

void OptFunctionSLRA::computePhat( gsl_vector* p, const gsl_vector* x ) {
  computeR(x, myTmpR);
  gsl_vector_memcpy(p, myFun.getP());
  myFun.computeCorrection(p, myTmpR);
}


