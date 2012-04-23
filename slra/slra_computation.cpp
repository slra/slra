#include <memory.h>

extern "C" {

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

}

#include "slra.h"

CostFunction::CostFunction( Structure *s, 
    int r, const gsl_vector *p, opt_and_info *opt, gsl_matrix *perm  ) : 
    myRank(r), myP(p), isGCD(opt->gcd) {
  myStruct = s; //.clone();     

  gsl_blas_ddot(myP, myP, &myPNorm);

  myGam = myStruct->createCholesky(getNplusD() - r, opt->reggamma);
  myDeriv = myStruct->createDGamma(getNplusD() - r);
      
  myMatr = gsl_matrix_alloc(getM(), getNplusD());
  myMatrMulPerm = gsl_matrix_alloc(getM(), getNplusD());
  myPerm = gsl_matrix_alloc(getNplusD(), getNplusD());
    
  myTmpGradR = gsl_matrix_alloc(getNplusD(), getD());
  myTmpGradR2 = gsl_matrix_alloc(getNplusD(), getD());

  myTmpR = gsl_matrix_alloc(getNplusD(), getD());
  myTmpYr = gsl_vector_alloc(getM() * getD());
  
  myTmpJacobianCol = gsl_vector_alloc(getM() * getD());

  myTmpGrad = gsl_matrix_alloc(myRank, getD());
  
  if (myStruct->getNp() != p->size) {
    throw new Exception("Inconsistent parameter vector\n");
  }

  myTmpCorr = gsl_vector_alloc(getNp());
  
  if (getM() < getNplusD()) {
    throw new Exception("Number of rows %d is less than the number of columns %d: ",
              getM(), getNplusD());
  }
  if (myStruct->getNp() < getM() * getD()) {
    throw new Exception("The inner minimization problem is overdetermined: " 
                        "m * (n-r) = %d, n_p = %d.\n", getM() * getD(), myStruct->getNp());
  }
    
  if (perm == NULL) {
    gsl_matrix_set_identity(myPerm);
  } else {
    if (perm->size1 != getNplusD() || perm->size2 != getNplusD()) {
      throw new Exception("Incorrect sizes of permutation matrix.\n");   
    }
    gsl_matrix_memcpy(myPerm, perm);
  }
  
  myStruct->fillMatrixFromP(myMatr, p);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myMatr, myPerm, 0.0, myMatrMulPerm);
}
  
CostFunction::~CostFunction() {
  delete myGam;
  delete myDeriv;

  gsl_matrix_free(myMatr);
  gsl_matrix_free(myMatrMulPerm);
  gsl_matrix_free(myPerm);
  gsl_matrix_free(myTmpGradR);
  gsl_matrix_free(myTmpGradR2);
  gsl_matrix_free(myTmpR);
  gsl_vector_free(myTmpYr);

  gsl_vector_free(myTmpJacobianCol);
  gsl_vector_free(myTmpCorr);

  gsl_matrix_free(myTmpGrad);
}



void CostFunction::computeR(gsl_matrix_const_view x_mat, gsl_matrix *R ) { 
  gsl_vector_view diag;
  gsl_matrix_view submat;
  int n = x_mat.matrix.size1, d = x_mat.matrix.size2;

  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  submat = gsl_matrix_submatrix(R, 0, 0, n, d); /* select x in (1,1) */
  gsl_matrix_memcpy(&submat.matrix, &x_mat.matrix); /* assign x */
  submat = gsl_matrix_submatrix(R, n, 0, d, d); /* select -I in (1,1)*/
  gsl_matrix_set_all(&submat.matrix,0);
  diag   = gsl_matrix_diagonal(&submat.matrix);     /* assign -I */
  gsl_vector_set_all(&diag.vector, -1);
  
  gsl_matrix_memcpy(myTmpGradR, R);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myPerm, myTmpGradR, 0.0, R);
}

void CostFunction::computeR( const gsl_vector * x, gsl_matrix *R ) { 
  computeR(gsl_matrix_const_view_vector(x, myRank, getD()), R);
}


void CostFunction::computeSr( gsl_matrix *R, gsl_vector *Sr ) {
  gsl_matrix_view SrMat = gsl_matrix_view_vector(Sr, getM(), getD()); 

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myMatr, R, 0.0, &SrMat.matrix);
}


void CostFunction::computeFuncAndGrad( const gsl_vector* x, double * f, gsl_vector *grad ) {
  computeRGammaSr(x, myTmpR, myTmpYr);

  if (f != NULL) {
    myGam->multiplyInvCholeskyVector(myTmpYr, 1);
    gsl_blas_ddot(myTmpYr, myTmpYr, f);
    
    if (isGCD) {
      *f = myPNorm - (*f);
    }
  }
  if (grad != NULL) {
    if (f != NULL) {
      myGam->multiplyInvCholeskyVector(myTmpYr, 0);
    } else {
      myGam->multiplyInvGammaVector(myTmpYr);
    }
    computeGradFromYr(myTmpYr, myTmpR, grad);
    if (isGCD) {
      gsl_vector_scale(grad, -1);
    }
  }
}

void CostFunction::computeFuncAndPseudoJacobianLs( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac ) {
  computeRGammaSr(x, myTmpR, myTmpYr);
  if (res != NULL) {
    myGam->multiplyInvCholeskyVector(myTmpYr, 1);
    gsl_vector_memcpy(res, myTmpYr);
  }
  if (jac != NULL) {  
    if (res != NULL) {
      myGam->multiplyInvCholeskyVector(myTmpYr, 0);
    } else {
      myGam->multiplyInvGammaVector(myTmpYr);      
    }
    computePseudoJacobianLsFromYr(myTmpYr, myTmpR, jac);
  } 
}


void CostFunction::computeCorrectionAndJacobian( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac ) {
  computeRGammaSr(x, myTmpR, myTmpYr);
  myGam->multiplyInvGammaVector(myTmpYr);

  if (res != NULL) {
    if (isGCD) {
      gsl_vector_memcpy(res, myP);
    } else {
      gsl_vector_set_zero(res);
    }
    myStruct->correctP(res, myTmpR, myTmpYr);
  }
  if (jac != NULL) {  
    computePseudoJacobianCorrectFromYr(myTmpYr, myTmpR, jac);
  } 
}


void CostFunction::computeJacobianZij( gsl_vector *res, int i, int j,
                                           gsl_vector* yr, gsl_matrix *R, double factor ) {
  myDeriv->calcDijGammaYr(res, R, myPerm, i, j, yr);
  gsl_vector_scale(res, -factor);
  
  for (int k = 0; k < getM(); k++) {  /* Convert to vector strides */
    (*gsl_vector_ptr(res, j + k * getD())) += gsl_matrix_get(myMatrMulPerm, k, i);
  }  
}

void CostFunction::computePseudoJacobianLsFromYr(  gsl_vector* yr, gsl_matrix *R, gsl_matrix *jac ) {
  int i, j;
//  gsl_vector m_col;

  for (i = 0; i < getN(); i++) {
    for (j = 0; j < getD(); j++) {
      computeJacobianZij(myTmpJacobianCol, i, j, yr, R);

      myGam->multiplyInvCholeskyVector(myTmpJacobianCol, 1);
//      m_col = gsl_matrix_column(jac, i * getD() + j).vector;
//      gsl_blas_dcopy(myTmpJacobianCol, &m_col);
      gsl_matrix_set_col(jac, i * getD() + j, myTmpJacobianCol);  
    }
  }
}

void CostFunction::computePseudoJacobianCorrectFromYr( gsl_vector* yr, gsl_matrix *R, gsl_matrix *jac ) {
  int i, j;
  gsl_vector_view jac_col, tmp_col;

  gsl_matrix_set_zero(jac);
  gsl_matrix_set_zero(myTmpGradR);

  for (i = 0; i < getN(); i++) {
    for (j = 0; j < getD(); j++) {
      jac_col = gsl_matrix_column(jac, i * getD() + j);

      gsl_vector_set_zero(myTmpCorr);

      /* Compute first term (correction of Gam^{-1} z_{ij}) */
      computeJacobianZij(myTmpJacobianCol, i, j, yr, R, 1);
      myGam->multiplyInvGammaVector(myTmpJacobianCol);
      myStruct->correctP(myTmpCorr, R, myTmpJacobianCol);

      /* Compute second term (gamma * dG_{ij} * yr) */ 
      gsl_matrix_set_zero(myTmpGradR);
      tmp_col = gsl_matrix_column(myPerm, i);
      gsl_matrix_set_col(myTmpGradR, j, &tmp_col.vector);
      myStruct->correctP(myTmpCorr, myTmpGradR, yr);

      /* Set to zero used column * /
      tmp_col = gsl_matrix_column(myTmpGradR, j);
      gsl_vector_set_zero(&tmp_col.vector); */
      
      gsl_vector_memcpy(&jac_col.vector, myTmpCorr);
    }
  }
}

void CostFunction::computeGradFromYr( gsl_vector* yr, gsl_matrix *R, gsl_vector *grad ) {
  gsl_matrix_view yr_matr = gsl_matrix_view_vector(yr, getM(), getD());
  gsl_matrix_view grad_matr = gsl_matrix_view_vector(grad, getN(), getD());
  gsl_matrix_view perm_sub_matr = gsl_matrix_submatrix(myPerm, 0, 0, getNplusD(), getN());

  /* Compute gradient of f(R) */ 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myMatr, &yr_matr.matrix, 0.0, myTmpGradR);
  myDeriv->calcYrtDgammaYr(myTmpGradR2, R, yr);
  gsl_matrix_sub(myTmpGradR, myTmpGradR2);

  /* Compute gradient of f_{\Phi}(X) */ 
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &perm_sub_matr.matrix, myTmpGradR, 0.0, &grad_matr.matrix);
}

void CostFunction::computeCorrection( gsl_vector* p, const gsl_vector* x ) {
  computeRGammaSr(x, myTmpR, myTmpYr);
  myGam->multiplyInvGammaVector(myTmpYr);
  myStruct->correctP(p, myTmpR, myTmpYr);
}

/*


class FnComputation {


  void computeLsPseudoJacobianFromYr( gsl_matrix *lsJacobian, const gsl_matrix *R, 
  {
    computePseudoJacobianFirstTerm(lsJacobian);
    computePseudoJacobianSecondTerm(P->tmpJacobian2, R, Phi);
    
    gsl_matrix_scale(P->tmpJacobian2, 0.5);
    gsl_matrix_sub(lsJacobian, P->tmpJacobian2);
  }


  void computeLsFun( gsl_vector *lsFun, const gsl_matrix *R ) {
    computeCholeskyOfGamma(R);
    computeYr(LsFunValue, R);
    multiplyInvCholesky(LsFunValue, 1);
  }


  void computeLsPseudoJacobianOnly( gsl_matrix *lsJacobian, const gsl_matrix *R, 
           const gsl_matrix *Phi ) {
    computeCholeskyOfGamma(R);
    computeYr(P->tmpYr, R);

    computeLsPseudoJacobianFromYr(computeLsPseudoJacobianFromYr, R, Phi);
  }

  void computeLsPseudoJacobianAndFun( gsl_vector *lsFun, gsl_matrix *lsJacobian, 
           const gsl_matrix *R, const gsl_matrix *Phi ) {
    computeCholeskyOfGamma(R);
    computeYr(P->tmpYr, R);
    
       
  }
}*/

