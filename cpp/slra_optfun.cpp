#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

OptFunctionSLRA::OptFunctionSLRA( CostFunction &fun, gsl_matrix *Phi, gsl_matrix *Psi ) : 
    myFun(fun) {
  if (Psi == NULL) {
    myPsi = NULL;
  } else {
    myPsi = gsl_matrix_alloc(Psi->size1, Psi->size2);
    gsl_matrix_memcpy(myPsi, Psi); 
  }
  if (Phi == NULL) {
    myPhi = NULL;
  } else {
    myPhi = gsl_matrix_alloc(Phi->size1, Phi->size2);
    gsl_matrix_memcpy(myPhi, Phi); 
  }
  
  
  if (Phi == NULL && Psi == NULL) {
    myPerm = gsl_matrix_alloc(myFun.getSt()->getM(), myFun.getSt()->getM());
    gsl_matrix_set_identity(myPerm);
  } else {
    if (Psi != NULL) {
      if (Phi != NULL) {
        myPerm = gsl_matrix_alloc(Phi->size1, Psi->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Phi, Psi, 0, myPerm); 
      } else {
        myPerm = gsl_matrix_alloc(Psi->size1, Psi->size2);
        gsl_matrix_memcpy(myPerm, Psi); 
      }
    } else {
      myPerm = gsl_matrix_alloc(Phi->size1, Phi->size2);
      gsl_matrix_memcpy(myPerm, Phi); 
    }
    if (myPerm->size1 != myFun.getSt()->getM() || myPerm->size1 < myPerm->size2) {
      throw new Exception("Incorrect sizes of permutation matrix.\n");   
    }
  }
  myTmpR = gsl_matrix_alloc(myFun.getSt()->getM(), myFun.getD());
  myTmpXId = gsl_matrix_alloc(myPerm->size2, myFun.getD());
}

OptFunctionSLRA::~OptFunctionSLRA()  {
  gsl_matrix_free(myTmpR);
  gsl_matrix_free(myPerm);
  gsl_matrix_free(myTmpXId);
  if (myPhi != NULL) {
    gsl_matrix_free(myPhi);
  }
  if (myPsi != NULL) {
    gsl_matrix_free(myPsi);
  }
}

int OptFunctionSLRA::getNvar() { 
  return getRank() * myFun.getD(); 
}

void OptFunctionSLRA::computeR( const gsl_vector * x, gsl_matrix *R ) { 
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector(x, getRank(),
                                    myFun.getD());

  X2XId(&x_mat.matrix, myTmpXId);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPerm, myTmpXId, 0, R);
}

void OptFunctionSLRA::computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad ) {
  computeR(x, myTmpR);

  if (grad == NULL) {
    myFun.computeFuncAndGrad(myTmpR, f, NULL, NULL);
  } else {
    gsl_matrix grad_matr = gsl_matrix_view_vector(grad, getRank(), myFun.getD()).matrix;
    gsl_matrix perm_sub_matr = gsl_matrix_submatrix(myPerm, 0, 0, 
                                      myFun.getSt()->getM(), getRank()).matrix;

    myFun.computeFuncAndGrad(myTmpR, f, myPerm, &grad_matr);
  }                   
}   

void OptFunctionSLRACholesky::computeFuncAndJac( const gsl_vector* x, 
         gsl_vector *res, gsl_matrix *jac ) {
  gsl_matrix perm_sub_matr = gsl_matrix_submatrix(myPerm, 0, 0, 
                                 myFun.getSt()->getM(), getRank()).matrix;         
  computeR(x, myTmpR);
  myFun.computeFuncAndPseudoJacobianLs(myTmpR, &perm_sub_matr, res, jac); 
}   

void OptFunctionSLRACorrection::computeFuncAndJac( const gsl_vector* x, 
         gsl_vector *res, gsl_matrix *jac ) {
  gsl_matrix perm_sub_matr = gsl_matrix_submatrix(myPerm, 0, 0, 
                                 myFun.getSt()->getM(), getRank()).matrix;         
  computeR(x, myTmpR);
  myFun.computeCorrectionAndJacobian(myTmpR, &perm_sub_matr, res, jac); 
}   


void OptFunctionSLRA::computeDefaultRTheta( gsl_matrix *RTheta ) {
  int c_size1 = myFun.getSt()->getN(), 
      c_size2 = (myPhi == NULL ? myFun.getSt()->getM(): myPhi->size2);
  size_t status = 0;
  size_t minus1 = -1;
  double tmp;

  gsl_matrix * tempc = gsl_matrix_alloc(c_size1, c_size2);
  if (myPhi == NULL) {
    gsl_matrix_memcpy(tempc, myFun.getSMatr());
  } else {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myFun.getSMatr(), myPhi, 0, tempc);
  }
  
  gsl_matrix * tempu = gsl_matrix_alloc(c_size2, c_size2);
  double *s = new double[mymin(c_size1, c_size2)];
  
  /* Determine optimal work */
  size_t lwork;
  dgesvd_("A", "N", &tempc->size2, &tempc->size1, tempc->data, &tempc->tda, s,
     tempu->data, &tempu->size2, NULL, &tempc->size1, &tmp, &minus1, &status);
  double *work = new double[(lwork = tmp)];
  /* Compute low-rank approximation */ 
  dgesvd_("A", "N", &tempc->size2, &tempc->size1, tempc->data, &tempc->tda, s,
     tempu->data, &tempu->size2, NULL, &tempc->size1, work, &lwork, &status);

  if (status) {
    delete [] s;  
    delete [] work;  
    gsl_matrix_free(tempc);
    gsl_matrix_free(tempu);
    throw new Exception("Error computing initial approximation: "
                        "DGESVD didn't converge\n");
  }

  gsl_matrix_transpose(tempu);
  gsl_matrix_view RlraT;
  RlraT = gsl_matrix_submatrix(tempu, 0, tempu->size2 - RTheta->size2, 
                               tempu->size1, RTheta->size2);
  gsl_matrix_memcpy(RTheta, &(RlraT.matrix));
    
  delete [] s;  
  delete [] work;  
  gsl_matrix_free(tempc);
  gsl_matrix_free(tempu);
}

void OptFunctionSLRA::X2XId( const gsl_matrix *x, gsl_matrix *XId ) { 
  gsl_vector diag;
  gsl_matrix sm;
  int n = x->size1, d = x->size2;

  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  gsl_matrix_memcpy(&(sm = gsl_matrix_submatrix(XId, 0, 0, n, d).matrix), x); 
  gsl_matrix_set_all(&(sm = gsl_matrix_submatrix(XId, n, 0, d, d).matrix), 0);
  gsl_vector_set_all(&(diag = gsl_matrix_diagonal(&sm).vector), -1);
}

void OptFunctionSLRA::PQ2XId( const gsl_matrix *R, gsl_matrix * x ) {
  size_t status = 0;
  gsl_matrix *tempR = gsl_matrix_alloc(R->size1, R->size2);
  gsl_matrix_memcpy(tempR, R);
  gsl_matrix_view RQt = gsl_matrix_submatrix(tempR, 0, 0, tempR->size1,
                                                          tempR->size2);
  /* Solve AX = B, where  R_theta = [B A] */
  gsl_matrix_view B =  gsl_matrix_submatrix(&(RQt.matrix), 0, 0, 
        RQt.matrix.size1 - RQt.matrix.size2, RQt.matrix.size2);
  gsl_matrix_view A =  gsl_matrix_submatrix(&(RQt.matrix), 
        RQt.matrix.size1 - RQt.matrix.size2, 0, 
        RQt.matrix.size2, RQt.matrix.size2);

  size_t *pivot = new size_t[A.matrix.size2];
  /* TODO: Check for singularity of A */
  dgesv_(&(A.matrix.size2), &(B.matrix.size1), A.matrix.data, &(A.matrix.tda), 
         pivot,  B.matrix.data, &(B.matrix.tda), &status);  
  delete [] pivot;
         
  gsl_matrix_memcpy(x, &(B.matrix));
  gsl_matrix_scale(x, -1.0);
  gsl_matrix_free(tempR);
}

void OptFunctionSLRA::RTheta2x( gsl_matrix *RTheta, gsl_vector *x ) {
  gsl_matrix x_mat = gsl_matrix_const_view_vector(x, 
        myTmpXId->size1 - myTmpXId->size2, myTmpXId->size2).matrix;
  if (myPsi == NULL) {
    gsl_matrix_memcpy(myTmpXId, RTheta);
  } else {
    ls_solve(myPsi, RTheta, myTmpXId);
  }

  PQ2XId(myTmpXId, &x_mat);
}

void OptFunctionSLRA::x2RTheta( gsl_matrix *RTheta, const gsl_vector *x ) {
  gsl_matrix x_mat = gsl_matrix_const_view_vector(x, 
        myTmpXId->size1 - myTmpXId->size2, myTmpXId->size2).matrix;
  X2XId(&x_mat, myTmpXId);

  if (myPsi == NULL) {
    gsl_matrix_memcpy(RTheta, myTmpXId);
  } else {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPsi, myTmpXId, 0, RTheta);
  }
}

void OptFunctionSLRA::computeDefaultx( gsl_vector *x ) {
 gsl_matrix *Rtheta;
 Rtheta = gsl_matrix_alloc(myPhi == NULL ? myFun.getSt()->getM() : 
                                           myPhi->size2 , myFun.getD());
 computeDefaultRTheta(Rtheta);
 RTheta2x(Rtheta, x);
 gsl_matrix_free(Rtheta);
}

void OptFunctionSLRA::computeCorrection( gsl_vector* p, const gsl_vector* x ) {
  computeR(x, myTmpR);
  myFun.computeCorrection(p, myTmpR);
}


