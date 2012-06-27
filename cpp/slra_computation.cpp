#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

CostFunction::CostFunction( Structure *s, int d, const gsl_vector *p, 
    OptimizationOptions *opt, gsl_matrix *perm ) :  
                                            myP(p), isGCD(opt->gcd) {
  myStruct = s;     

  gsl_blas_ddot(myP, myP, &myPNorm);

  if (perm == NULL) {
    myPerm = gsl_matrix_alloc(getM(), getM());
  } else {
    myPerm = gsl_matrix_alloc(perm->size1, perm->size2);
  }
  
  myRank = myPerm->size2 - d;

  myGam = myStruct->createCholesky(getD(), opt->reggamma);
  myDeriv = myStruct->createDGamma(getD());

  myMatr = gsl_matrix_alloc(getN(), getM());
  myMatrMulPerm = gsl_matrix_alloc(getN(), myPerm->size2);
  myTmpThetaExt = gsl_matrix_alloc(myPerm->size2, getD());

  myTmpGradR = gsl_matrix_alloc(getM(), getD());
  myTmpGradR2 = gsl_matrix_alloc(getM(), getD());

  myTmpR = gsl_matrix_alloc(getM(), getD());
  myTmpYr = gsl_vector_alloc(getN() * getD());
  myTmpJacobianCol = gsl_vector_alloc(getN() * getD());

  myTmpGrad = gsl_matrix_alloc(myRank, getD());

  myEye = gsl_matrix_alloc(getM(), getM());
  gsl_matrix_set_identity(myEye);
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
    
  if (perm == NULL) {
    gsl_matrix_set_identity(myPerm);
  } else {
    if (perm->size1 != getM() || perm->size1 < perm->size2) {
      throw new Exception("Incorrect sizes of permutation matrix.\n");   
    }
    gsl_matrix_memcpy(myPerm, perm);
  }
 
  myStruct->fillMatrixFromP(myMatr, p);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myMatr, myPerm,
                 0.0, myMatrMulPerm);
}
  
CostFunction::~CostFunction() {
  delete myGam;
  delete myDeriv;

  gsl_matrix_free(myMatr);
  gsl_matrix_free(myMatrMulPerm);
  gsl_matrix_free(myPerm);
  gsl_matrix_free(myTmpThetaExt);
  gsl_matrix_free(myTmpGradR);
  gsl_matrix_free(myTmpGradR2);
  gsl_matrix_free(myTmpR);
  
  gsl_matrix_free(myEye);
  gsl_matrix_free(myTmpJac);

  gsl_vector_free(myTmpYr);

  gsl_vector_free(myTmpJacobianCol);
  gsl_vector_free(myTmpCorr);

  gsl_matrix_free(myTmpGrad);
}

void CostFunction::computeRGammaSr( const gsl_vector *x, gsl_matrix *R, 
         gsl_vector *Sr, bool regularize_gamma ) {
  computeR(x, myTmpR);
  myGam->calcGammaCholesky(myTmpR, regularize_gamma);
  computeSr(myTmpR, Sr);
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

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myPerm, myTmpThetaExt, 0.0,
      R);
}

void CostFunction::computeR( const gsl_vector * x, gsl_matrix *R ) { 
  computeR(gsl_matrix_const_view_vector(x, getRank(), getD()), R);
}

void CostFunction::computeSr( gsl_matrix *R, gsl_vector *Sr ) {
  gsl_matrix_view SrMat = gsl_matrix_view_vector(Sr, getN(), getD()); 

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myMatr, R, 0.0, 
                 &SrMat.matrix);
}

void CostFunction::computeFuncAndGrad( const gsl_vector* x, double * f, 
                                       gsl_vector *grad ) {
  computeRGammaSr(x, myTmpR, myTmpYr);

  if (f != NULL) {
    myGam->multInvCholeskyVector(myTmpYr, 1);
    gsl_blas_ddot(myTmpYr, myTmpYr, f);
    if (isGCD) {
      *f = myPNorm - (*f);
    }
  }
  if (grad != NULL) {
    if (f != NULL) {
      myGam->multInvCholeskyVector(myTmpYr, 0);
    } else {
      myGam->multInvGammaVector(myTmpYr);
    }
    computeGradFromYr(myTmpYr, myTmpR, grad);
    if (isGCD) {
      gsl_vector_scale(grad, -1);
    }
  }
}

void CostFunction::computeFuncAndPseudoJacobianLs( const gsl_vector* x, 
         gsl_vector *res, gsl_matrix *jac ) {
  computeRGammaSr(x, myTmpR, myTmpYr);

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
  computeRGammaSr(x, myTmpR, myTmpYr);
  myGam->multInvGammaVector(myTmpYr);

  if (res != NULL) {
    if (isGCD) {
      gsl_vector_memcpy(res, myP);
    } else {
      gsl_vector_set_zero(res);
    }
    myStruct->correctP(res, myTmpR, myTmpYr, 1);
  }
  if (jac != NULL) {  
    computeJacobianOfCorrection(myTmpYr, myTmpR, jac);
  } 
}


void CostFunction::computeJacobianZijOld( gsl_vector *res, int i, int j,
         gsl_vector* yr, gsl_matrix *R, double factor ) {
  myDeriv->calcDijGammaYr(res, R, myPerm, i, j, yr);
  gsl_vector_scale(res, -factor);
  
  for (int k = 0; k < getN(); k++) {  /* Convert to vector strides */
    (*gsl_vector_ptr(res, j + k * getD())) += 
          gsl_matrix_get(myMatrMulPerm, k, i);
  }  
}


void CostFunction::computeJacobianZij( gsl_vector *res, int i, int j,
         gsl_vector* yr, gsl_matrix *R, double factor ) {
  myDeriv->calcDijGammaYr(res, R, myEye, i, j, yr);
  gsl_vector_scale(res, -factor);
  
  for (int k = 0; k < getN(); k++) {  /* Convert to vector strides */
    (*gsl_vector_ptr(res, j + k * getD())) += 
          gsl_matrix_get(myMatr, k, i);
  }  
}


void CostFunction::computePseudoJacobianLsFromYrOld( gsl_vector* yr, 
         gsl_matrix *R, gsl_matrix *jac ) {
  for (size_t i = 0; i < getRank(); i++) {
    for (size_t j = 0; j < getD(); j++) {
      computeJacobianZij(myTmpJacobianCol, i, j, yr, R);

      myGam->multInvCholeskyVector(myTmpJacobianCol, 1);
      gsl_matrix_set_col(jac, i * getD() + j, myTmpJacobianCol);  
    }
  }
}

void CostFunction::computePseudoJacobianLsFromYr( gsl_vector* yr, 
         gsl_matrix *R, gsl_matrix *jac ) {
  size_t i, j;
  for (i = 0; i < getM(); i++) {
    for (j = 0; j < getD(); j++) {
      gsl_vector tmpJacRow = gsl_matrix_row(myTmpJac, i + j * getM()).vector;
    
      computeJacobianZij(&tmpJacRow, i, j, yr, R);
    }
  }
  
  for (i = 0; i < getRank(); i++) {
    for (j = 0; j < getD(); j++) {
      gsl_matrix subJ = gsl_matrix_submatrix(myTmpJac, j * getM(), 0, getM(),
                            myTmpJac->size2).matrix;
      gsl_vector phiCol = gsl_matrix_column(myPerm, i).vector;
      gsl_blas_dgemv(CblasTrans, 1.0, &subJ, &phiCol, 0.0, myTmpJacobianCol);
    
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

  for (i = 0; i < getRank(); i++) {
    for (j = 0; j < getD(); j++) {
      jac_col = gsl_matrix_column(jac, i * getD() + j);
      gsl_vector_set_zero(myTmpCorr);

      /* Compute first term (correction of Gam^{-1} z_{ij}) */
      computeJacobianZij(myTmpJacobianCol, i, j, yr, R, 1);
      myGam->multInvGammaVector(myTmpJacobianCol);
      myStruct->correctP(myTmpCorr, R, myTmpJacobianCol, 1);

      /* Compute second term (gamma * dG_{ij} * yr) */ 
      gsl_matrix_set_zero(myTmpGradR);
      tmp_col = gsl_matrix_column(myPerm, i);
      gsl_matrix_set_col(myTmpGradR, j, &tmp_col.vector);
      myStruct->correctP(myTmpCorr, myTmpGradR, yr, 1);

      /* Set to zero used column * /
      tmp_col = gsl_matrix_column(myTmpGradR, j);
      gsl_vector_set_zero(&tmp_col.vector); */
      
      gsl_vector_memcpy(&jac_col.vector, myTmpCorr);
    }
  }
}

void CostFunction::computeGradFromYr( gsl_vector* yr, gsl_matrix *R, 
                                      gsl_vector *grad ) {
  gsl_matrix_view yr_matr = gsl_matrix_view_vector(yr, getN(), getD());
  gsl_matrix_view grad_matr = gsl_matrix_view_vector(grad, getRank(), getD());
  gsl_matrix_view perm_sub_matr = gsl_matrix_submatrix(myPerm, 0, 0, 
                                      getM(), getRank());

  /* Compute gradient of f(R) */ 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myMatr, &yr_matr.matrix,
                 0.0, myTmpGradR);
  myDeriv->calcYrtDgammaYr(myTmpGradR2, R, yr);
  gsl_matrix_sub(myTmpGradR, myTmpGradR2);

  /* Compute gradient of f_{\Phi}(X) */ 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &perm_sub_matr.matrix, 
                 myTmpGradR, 0.0, &grad_matr.matrix);
}

void CostFunction::computeCorrection( gsl_vector* p, const gsl_vector* x ) {
  try  {
    computeRGammaSr(x, myTmpR, myTmpYr, false);
  } catch (Exception *e) {
    if (!strncmp(e->getMessage(), "Gamma", 5)) {
      delete e;
      e = new Exception("Gamma matrix is singular. "
                        "Unable to compute the correction.\n");
    }
    throw e;
  }
  myGam->multInvGammaVector(myTmpYr);
  if (isGCD) {
    gsl_vector_set_zero(p);
    myStruct->correctP(p, myTmpR, myTmpYr, 0);
    gsl_vector_scale(p, -1.0);
  } else {
    myStruct->correctP(p, myTmpR, myTmpYr);
  }
}

void CostFunction::computeDefaultRTheta( gsl_matrix *RTheta ) {
  const gsl_matrix *c = getPhiSMatr();
  size_t status = 0;
  size_t minus1 = -1;
  double tmp;

  gsl_matrix * tempc = gsl_matrix_alloc(c->size1, c->size2);
  gsl_matrix_memcpy(tempc, c);
  gsl_matrix * tempu = gsl_matrix_alloc(c->size2, c->size2);
  double *s = new double[mymin(c->size1, c->size2)];
  
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



