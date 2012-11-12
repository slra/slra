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

CostFunction::CostFunction( Structure *s, int d, const gsl_vector *p, 
                            gsl_matrix *perm, double reggamma ) :  
                            myP(p), myD(d), myStruct(s) {
  myGam = myStruct->createCholesky(getD(), reggamma);
  myDeriv = myStruct->createDGamma(getD());

  myMatr = gsl_matrix_alloc(getN(), getM());
  myTmpGradR = gsl_matrix_alloc(getM(), getD());

  myTmpR = gsl_matrix_alloc(getM(), getD());
  myTmpYr = gsl_vector_alloc(getN() * getD());
  myTmpJacobianCol = gsl_vector_alloc(getN() * getD());
  myTmpJac = gsl_matrix_alloc(getM() * getD(), getN() * getD());

  if (myStruct->getNp() > p->size) {
    throw new Exception("Inconsistent parameter vector\n");
  }
  
  myTmpCorr = gsl_vector_alloc(getNp());
  
  if (getN() < getM()) {
    throw new Exception("Number of columns %d is less than "
                        "the number of rows %d.", getN(), getM());
  }

  if (myStruct->getNp() < getN() * getD()) {
    throw new Exception("The inner minimization problem is overdetermined: " 
        "n * (m-r) = %d, n_p = %d.\n", getN() * getD(), myStruct->getNp());
  }
  myStruct->fillMatrixFromP(myMatr, p);
    
  if (perm == NULL) {
    myPerm = gsl_matrix_alloc(getM(), getM());
    gsl_matrix_set_identity(myPerm);
  } else {
    if (perm->size1 != getM() || perm->size1 < perm->size2) {
      throw new Exception("Incorrect sizes of permutation matrix.\n");   
    }
    myPerm = gsl_matrix_alloc(perm->size1, perm->size2);
    gsl_matrix_memcpy(myPerm, perm);
  }
  myTmpThetaExt = gsl_matrix_alloc(myPerm->size2, getD());
}
  
CostFunction::~CostFunction() {
  delete myGam;
  delete myDeriv;

  gsl_matrix_free(myMatr);
  gsl_matrix_free(myTmpGradR);
  gsl_matrix_free(myTmpR);
  gsl_vector_free(myTmpYr);

  gsl_matrix_free(myTmpJac);
  gsl_vector_free(myTmpJacobianCol);
  gsl_vector_free(myTmpCorr);

  gsl_matrix_free(myPerm);
  gsl_matrix_free(myTmpThetaExt);
}


void CostFunction::X2Rtheta( const gsl_matrix *x, gsl_matrix *RTheta ) { 
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

void CostFunction::computeR(gsl_matrix_const_view x_mat, gsl_matrix *R ) { 
  X2Rtheta(&x_mat.matrix, myTmpThetaExt);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPerm, myTmpThetaExt, 0, R);
}

void CostFunction::computeR( const gsl_vector * x, gsl_matrix *R ) { 
  computeR(gsl_matrix_const_view_vector(x, getRank(), getD()), R);
}

void CostFunction::computeSr( const gsl_matrix *R, gsl_vector *Sr ) {
  gsl_matrix_view SrMat = gsl_matrix_view_vector(Sr, getN(), getD()); 
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myMatr, R, 0, &SrMat.matrix);
}

void CostFunction::computeGammaSr( const gsl_matrix *R, gsl_vector *Sr, bool regularize_gamma ) {
  myGam->calcGammaCholesky(R, regularize_gamma);
  computeSr(R, Sr);
} 

void CostFunction::computeFuncAndGrad( const gsl_matrix* R, double * f, 
                                       gsl_matrix *gradR ) {
  computeGammaSr(R, myTmpYr);

  if (f != NULL) {
    myGam->multInvCholeskyVector(myTmpYr, 1);
    gsl_blas_ddot(myTmpYr, myTmpYr, f);
  }
  if (gradR != NULL) {
    if (f != NULL) {
      myGam->multInvCholeskyVector(myTmpYr, 0);
    } else {
      myGam->multInvGammaVector(myTmpYr);
    }
    gsl_matrix_view yr_matr = gsl_matrix_view_vector(myTmpYr, getN(), getD());
    myDeriv->calcYrtDgammaYr(gradR, R, myTmpYr);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myMatr, &yr_matr.matrix,
                   -1.0, gradR);
  }
}

void CostFunction::computeFuncAndPseudoJacobianLs( const gsl_vector* x, 
         gsl_vector *res, gsl_matrix *jac ) {
  computeR(x, myTmpR);
  computeGammaSr(myTmpR, myTmpYr);

  if (res != NULL) {
    myGam->multInvCholeskyVector(myTmpYr, 1);
    gsl_vector_memcpy(res, myTmpYr);
  }

  if (jac != NULL) {  
    if (res != NULL) {
      myGam->multInvCholeskyVector(myTmpYr, 0);
    } else {
      myGam->multInvGammaVector(myTmpYr);      
    }
    computePseudoJacobianLsFromYr(myTmpYr, myTmpR, jac);
  } 
}

void CostFunction::computeCorrectionAndJacobian( const gsl_vector* x, 
         gsl_vector *res, gsl_matrix *jac ) {
  computeR(x, myTmpR);
  computeGammaSr(myTmpR, myTmpYr);
  myGam->multInvGammaVector(myTmpYr);

  if (res != NULL) {
    gsl_vector_set_zero(res);
    myStruct->correctP(res, myTmpR, myTmpYr, 1);
  }
  if (jac != NULL) {  
    computeJacobianOfCorrection(myTmpYr, myTmpR, jac);
  } 
}

void CostFunction::computeZmatTmpJac( gsl_vector* yr, gsl_matrix *R, double factor ) {
  size_t i, j, k;
  for (i = 0; i < getM(); i++) {
    for (j = 0; j < getD(); j++) {
      gsl_vector tmpJacRow = gsl_matrix_row(myTmpJac, i + j * getM()).vector;
    
      myDeriv->calcDijGammaYr(&tmpJacRow, R, i, j, yr);
      gsl_vector_scale(&tmpJacRow, -factor);
  
      for (int k = 0; k < getN(); k++) {  /* Convert to vector strides */
        (*gsl_vector_ptr(&tmpJacRow, j + k * getD())) += 
             gsl_matrix_get(myMatr, k, i);
      }  
    }
  }
}

void CostFunction::mulZmatPerm( gsl_vector* res, gsl_matrix *perm, size_t i, size_t j ) {
  gsl_matrix subJ = gsl_matrix_submatrix(myTmpJac, j * getM(), 0, getM(),
                         myTmpJac->size2).matrix;
  gsl_vector phiCol = gsl_matrix_column(perm, i).vector;
  gsl_blas_dgemv(CblasTrans, 1.0, &subJ, &phiCol, 0.0, res);
}

void CostFunction::computePseudoJacobianLsFromYr( gsl_vector* yr, 
         gsl_matrix *R, gsl_matrix *jac ) {
  computeZmatTmpJac(yr, R);
  
  for (size_t i = 0; i < getRank(); i++) {
    for (size_t j = 0; j < getD(); j++) {
      mulZmatPerm(myTmpJacobianCol, myPerm, i, j);
      myGam->multInvCholeskyVector(myTmpJacobianCol, 1);
      gsl_matrix_set_col(jac, i * getD() + j, myTmpJacobianCol);  
    }
  }
}

void CostFunction::computeJacobianOfCorrection( gsl_vector* yr, 
         gsl_matrix *R, gsl_matrix *jac ) {
  int i, j;
  gsl_vector_view jac_col, tmp_col;

  gsl_matrix_set_zero(jac);
  gsl_matrix_set_zero(myTmpGradR);

  computeZmatTmpJac(yr, R);

  for (i = 0; i < getRank(); i++) {
    for (j = 0; j < getD(); j++) {
      jac_col = gsl_matrix_column(jac, i * getD() + j);
      gsl_vector_set_zero(myTmpCorr);

      /* Compute first term (correction of Gam^{-1} z_{ij}) */
      mulZmatPerm(myTmpJacobianCol, myPerm, i, j);
      myGam->multInvGammaVector(myTmpJacobianCol);
      myStruct->correctP(myTmpCorr, R, myTmpJacobianCol, 1);

      /* Compute second term (gamma * dG_{ij} * yr) */ 
      gsl_matrix_set_zero(myTmpGradR);
      tmp_col = gsl_matrix_column(myPerm, i);
      gsl_matrix_set_col(myTmpGradR, j, &tmp_col.vector);
      myStruct->correctP(myTmpCorr, myTmpGradR, yr, 1);

      gsl_vector_memcpy(&jac_col.vector, myTmpCorr);
    }
  }
}

void CostFunction::computeCorrection( gsl_vector* p, const gsl_vector* x ) {
  try  {
    computeR(x, myTmpR);
    computeGammaSr(myTmpR, myTmpYr, false);
  } catch (Exception *e) {
    if (!strncmp(e->getMessage(), "Gamma", 5)) {
      delete e;
      e = new Exception("Gamma matrix is singular. "
                        "Unable to compute the correction.\n");
    }
    throw e;
  }
  myGam->multInvGammaVector(myTmpYr);
  myStruct->correctP(p, myTmpR, myTmpYr);
}

void CostFunction::computeDefaultRTheta( gsl_matrix *RTheta ) {
  int c_size1 = getN(), c_size2 = getPerm()->size2;
  size_t status = 0;
  size_t minus1 = -1;
  double tmp;

  gsl_matrix * tempc = gsl_matrix_alloc(c_size1, c_size2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, getSMatr(), getPerm(), 0, tempc);
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

void CostFunction::Rtheta2X( const gsl_matrix *R, gsl_matrix * x ) {
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



