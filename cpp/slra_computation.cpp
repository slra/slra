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
                            double reggamma ) :  
                            myP(p), myD(d), myStruct(s) {
  myGam = myStruct->createCholesky(getD(), reggamma);
  myDeriv = myStruct->createDGamma(getD());

  myMatr = gsl_matrix_alloc(myStruct->getN(), myStruct->getM());
  myTmpGradR = gsl_matrix_alloc(myStruct->getM(), getD());

  myTmpYr = gsl_vector_alloc(myStruct->getN() * getD());
  myTmpJacobianCol = gsl_vector_alloc(myStruct->getN() * getD());
  myTmpJac = gsl_matrix_alloc(myStruct->getM() * getD(), myStruct->getN() * getD());

  if (myStruct->getNp() > p->size) {
    throw new Exception("Inconsistent parameter vector\n");
  }
  
  myTmpCorr = gsl_vector_alloc(myStruct->getNp());
  
  if (myStruct->getN() < myStruct->getM()) {
    throw new Exception("Number of columns %d is less than "
                        "the number of rows %d.", myStruct->getN(), myStruct->getM());
  }

  if (myStruct->getNp() < myStruct->getN() * getD()) {
    throw new Exception("The inner minimization problem is overdetermined: " 
        "n * (m-r) = %d, n_p = %d.\n", myStruct->getN() * getD(), myStruct->getNp());
  }
  myStruct->fillMatrixFromP(myMatr, p);
}
  
CostFunction::~CostFunction() {
  delete myGam;
  delete myDeriv;

  gsl_matrix_free(myMatr);
  gsl_matrix_free(myTmpGradR);
  gsl_vector_free(myTmpYr);

  gsl_matrix_free(myTmpJac);
  gsl_vector_free(myTmpJacobianCol);
  gsl_vector_free(myTmpCorr);
}

void CostFunction::computeGammaSr( const gsl_matrix *R, gsl_vector *Sr, bool regularize_gamma ) {
  myGam->calcGammaCholesky(R, regularize_gamma);
  gsl_matrix_view SrMat = gsl_matrix_view_vector(Sr, getSt()->getN(), getD()); 
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myMatr, R, 0, &SrMat.matrix);
} 

void CostFunction::computeFuncAndGrad( const gsl_matrix* R, double * f, 
                                       gsl_matrix *perm, gsl_matrix *gradR ) {
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
    gsl_matrix_view yr_matr = gsl_matrix_view_vector(myTmpYr, getSt()->getN(), getD());
    myDeriv->calcYrtDgammaYr(myTmpGradR, R, myTmpYr);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myMatr, &yr_matr.matrix,
                   -1.0, myTmpGradR);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, perm, myTmpGradR, 0.0, gradR);
  }
}

void CostFunction::computeFuncAndPseudoJacobianLs( gsl_matrix *R, gsl_matrix *perm,  
                                                   gsl_vector *res, gsl_matrix *jac ) {
  computeGammaSr(R, myTmpYr);                                                   

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
    computePseudoJacobianLsFromYr(myTmpYr, R, perm, jac);
  } 
}

void CostFunction::computeCorrectionAndJacobian( gsl_matrix *R, gsl_matrix *perm, 
                                                 gsl_vector *res, gsl_matrix *jac ) {
  computeGammaSr(R, myTmpYr);                                                   
  myGam->multInvGammaVector(myTmpYr);
  if (res != NULL) {
    gsl_vector_set_zero(res);
    myStruct->correctP(res, R, myTmpYr, 1);
  }
  if (jac != NULL) {  
    computeJacobianOfCorrection(myTmpYr, R, perm, jac);
  } 
}

void CostFunction::computeZmatTmpJac( gsl_vector* yr, gsl_matrix *R, double factor ) {
  size_t i, j, k;
  for (i = 0; i < getSt()->getM(); i++) {
    for (j = 0; j < getD(); j++) {
      gsl_vector tmpJacRow = gsl_matrix_row(myTmpJac, i + j * getSt()->getM()).vector;
    
      myDeriv->calcDijGammaYr(&tmpJacRow, R, i, j, yr);
      gsl_vector_scale(&tmpJacRow, -factor);
  
      for (int k = 0; k < getSt()->getN(); k++) {  /* Convert to vector strides */
        (*gsl_vector_ptr(&tmpJacRow, j + k * getD())) += 
             gsl_matrix_get(myMatr, k, i);
      }  
    }
  }
}

void CostFunction::mulZmatPerm( gsl_vector* res, gsl_matrix *perm, size_t i, size_t j ) {
  gsl_matrix subJ = gsl_matrix_submatrix(myTmpJac, j * getSt()->getM(), 0, getSt()->getM(),
                         myTmpJac->size2).matrix;
  gsl_vector phiCol = gsl_matrix_column(perm, i).vector;
  gsl_blas_dgemv(CblasTrans, 1.0, &subJ, &phiCol, 0.0, res);
}

void CostFunction::computePseudoJacobianLsFromYr( gsl_vector* yr, 
         gsl_matrix *R, gsl_matrix *perm, gsl_matrix *jac ) {
  computeZmatTmpJac(yr, R);
  
  for (size_t i = 0; i < perm->size2; i++) {
    for (size_t j = 0; j < getD(); j++) {
      mulZmatPerm(myTmpJacobianCol, perm, i, j);
      myGam->multInvCholeskyVector(myTmpJacobianCol, 1);
      gsl_matrix_set_col(jac, i * getD() + j, myTmpJacobianCol);  
    }
  }
}

void CostFunction::computeJacobianOfCorrection( gsl_vector* yr, 
         gsl_matrix *R, gsl_matrix *perm, gsl_matrix *jac ) {
  int i, j;
  gsl_vector_view jac_col, tmp_col;

  gsl_matrix_set_zero(jac);
  gsl_matrix_set_zero(myTmpGradR);

  computeZmatTmpJac(yr, R);

  for (i = 0; i < perm->size2; i++) {
    for (j = 0; j < getD(); j++) {
      jac_col = gsl_matrix_column(jac, i * getD() + j);
      gsl_vector_set_zero(myTmpCorr);

      /* Compute first term (correction of Gam^{-1} z_{ij}) */
      mulZmatPerm(myTmpJacobianCol, perm, i, j);
      myGam->multInvGammaVector(myTmpJacobianCol);
      myStruct->correctP(myTmpCorr, R, myTmpJacobianCol, 1);

      /* Compute second term (gamma * dG_{ij} * yr) */ 
      gsl_matrix_set_zero(myTmpGradR);
      tmp_col = gsl_matrix_column(perm, i);
      gsl_matrix_set_col(myTmpGradR, j, &tmp_col.vector);
      myStruct->correctP(myTmpCorr, myTmpGradR, yr, 1);

      gsl_vector_memcpy(&jac_col.vector, myTmpCorr);
    }
  }
}

void CostFunction::computeCorrection( gsl_vector* p, gsl_matrix *R ) {
  try  {
    computeGammaSr(R, myTmpYr, false);
  } catch (Exception *e) {
    if (!strncmp(e->getMessage(), "Gamma", 5)) {
      delete e;
      e = new Exception("Gamma matrix is singular. "
                        "Unable to compute the correction.\n");
    }
    throw e;
  }
  myGam->multInvGammaVector(myTmpYr);
  myStruct->correctP(p, R, myTmpYr);
}





