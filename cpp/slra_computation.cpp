#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

CostFunction::CostFunction( const gsl_vector *p, Structure *s, int d, 
                         gsl_matrix *Phi ) : myP(NULL), myD(d), myStruct(s), 
                         myReggamma(SLRA_DEF_reggamma) {
  if (myStruct->getNp() > p->size) {
    throw new Exception("Inconsistent parameter vector\n");
  }
  if (myStruct->getN() < myStruct->getM()) {
    throw new Exception("Number of columns %d is less than "
                        "the number of rows %d.", myStruct->getN(), myStruct->getM());
  }
  if (myStruct->getNp() < myStruct->getN() * getD()) {
    throw new Exception("The inner minimization problem is overdetermined: " 
        "n * (m-r) = %d, n_p = %d.\n", myStruct->getN() * getD(), myStruct->getNp());
  }

  if (Phi != NULL) {
    if (Phi->size1 != myStruct->getM() || Phi->size2 > Phi->size1) {
      throw new Exception("Incompatible Phi matrix\n");
    }
    myPhi = gsl_matrix_alloc(Phi->size1, Phi->size2);
    gsl_matrix_memcpy(myPhi, Phi);
  } else {
    myPhi = gsl_matrix_alloc(myStruct->getM(), myStruct->getM());
    gsl_matrix_set_identity(myPhi);
  }

  if (d >= getNrow() || d <= 0) {
    throw new Exception("Incorrect rank given\n");
  }

  myP = gsl_vector_alloc(p->size);
  gsl_vector_memcpy(myP, p);
  myRorig = gsl_matrix_alloc(getM(), getD());
  myPhiPermCol = gsl_vector_alloc(getM());
  myGam = myStruct->createCholesky(getD());
  myDeriv = myStruct->createDGamma(getD());
  myMatr = gsl_matrix_alloc(myStruct->getN(), myStruct->getM());
  myTmpGradR = gsl_matrix_alloc(getM(), getD());
  myTmpGradR2 = gsl_matrix_alloc(getNrow(), getD());
  myTmpYr = gsl_vector_alloc(myStruct->getN() * getD());
  myTmpJacobianCol = gsl_vector_alloc(myStruct->getN() * getD());
  myTmpJac = gsl_matrix_alloc(myStruct->getM() * getD(), myStruct->getN() * getD());
  myTmpCorr = gsl_vector_alloc(myStruct->getNp());
  myStruct->fillMatrixFromP(myMatr, p);
}
  
CostFunction::~CostFunction() {
  delete myGam;
  delete myDeriv;
  gsl_vector_free(myP);
  gsl_matrix_free(myPhi);
  gsl_matrix_free(myRorig);
  gsl_vector_free(myPhiPermCol);
  gsl_matrix_free(myMatr);
  gsl_matrix_free(myTmpGradR);
  gsl_matrix_free(myTmpGradR2);
  gsl_vector_free(myTmpYr);
  gsl_matrix_free(myTmpJac);
  gsl_vector_free(myTmpJacobianCol);
  gsl_vector_free(myTmpCorr);
}

void CostFunction::computeGammaSr( const gsl_matrix *R, 
                      gsl_matrix *Rorig, gsl_vector *Sr, bool regularize_gamma ) {
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPhi, R, 0, Rorig);
  myGam->calcGammaCholesky(Rorig, regularize_gamma ? myReggamma : 0);
  gsl_matrix_view SrMat = gsl_matrix_view_vector(Sr, getN(), getD()); 
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myMatr, Rorig, 0, &SrMat.matrix);
} 

void CostFunction::computeZmatTmpJac( gsl_vector* yr, gsl_matrix *Rorig, double factor ) {
  for (size_t i = 0; i < getM(); i++) {
    for (size_t j = 0; j < getD(); j++) {
      gsl_vector tJr = gsl_matrix_row(myTmpJac, i + j * getM()).vector;
    
      myDeriv->calcDijGammaYr(&tJr, Rorig, i, j, yr);
      gsl_vector_scale(&tJr, -factor);
      for (size_t k = 0; k < getN(); k++) {  /* Convert to vector strides */
        (*gsl_vector_ptr(&tJr, j + k * getD())) += gsl_matrix_get(myMatr, k, i);
      }  
    }
  }
}

void CostFunction::mulZmatPerm( gsl_vector* res, gsl_matrix *perm, size_t i, size_t j ) {
  gsl_matrix subJ = gsl_matrix_submatrix(myTmpJac, j * getM(), 0, getM(),
                         myTmpJac->size2).matrix;
  gsl_vector permCol = gsl_matrix_column(perm, i).vector;
  gsl_blas_dgemv(CblasNoTrans, 1.0, myPhi, &permCol, 0.0, myPhiPermCol);
  gsl_blas_dgemv(CblasTrans, 1.0, &subJ, myPhiPermCol, 0.0, res);
}

void CostFunction::computePseudoJacobianLsFromYr( gsl_vector* yr, 
         gsl_matrix *Rorig, gsl_matrix *perm, gsl_matrix *jac ) {
  computeZmatTmpJac(yr, Rorig);
  
  for (size_t i = 0; i < perm->size2; i++) {
    for (size_t j = 0; j < getD(); j++) {
      mulZmatPerm(myTmpJacobianCol, perm, i, j);
      myGam->multInvCholeskyVector(myTmpJacobianCol, 1);
      gsl_matrix_set_col(jac, i * getD() + j, myTmpJacobianCol);  
    }
  }
}

void CostFunction::computeJacobianOfCorrection( gsl_vector* yr, 
         gsl_matrix *Rorig, gsl_matrix *perm, gsl_matrix *jac ) {
  gsl_vector_view jac_col, tmp_col;
  gsl_matrix_set_zero(jac);
  gsl_matrix_set_zero(myTmpGradR);

  computeZmatTmpJac(yr, Rorig, 11);

  for (size_t i = 0; i < perm->size2; i++) {
    for (size_t j = 0; j < getD(); j++) {
      jac_col = gsl_matrix_column(jac, i * getD() + j);
      gsl_vector_set_zero(myTmpCorr);

      /* Compute first term (correction of Gam^{-1} z_{ij}) */
      mulZmatPerm(myTmpJacobianCol, perm, i, j);
      myGam->multInvGammaVector(myTmpJacobianCol);
      myStruct->correctP(myTmpCorr, Rorig, myTmpJacobianCol, 1);

      /* Compute second term (gamma * dG_{ij} * yr) */ 
      gsl_matrix_set_zero(myTmpGradR);
      gsl_vector permCol = gsl_matrix_column(perm, i).vector;
      gsl_blas_dgemv(CblasNoTrans, 1.0, myPhi, &permCol, 0.0, myPhiPermCol);
      gsl_matrix_set_col(myTmpGradR, j, myPhiPermCol);
      myStruct->correctP(myTmpCorr, myTmpGradR, yr, 1);

      gsl_vector_memcpy(&jac_col.vector, myTmpCorr);
    }
  }
}

void CostFunction::computeGradFromYr( gsl_vector* yr, const gsl_matrix *Rorig, 
                                      gsl_matrix *perm, gsl_matrix *gradR ) {
  gsl_matrix_view yr_matr = gsl_matrix_view_vector(yr, getN(), getD());
  myDeriv->calcYrtDgammaYr(myTmpGradR, Rorig, yr);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myMatr, &yr_matr.matrix,
                -1.0, myTmpGradR);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, myPhi, myTmpGradR,
                 0.0, myTmpGradR2);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, perm, myTmpGradR2, 0.0, gradR);
}

void CostFunction::computeFuncAndGrad( const gsl_matrix* R, double * f, 
                                       gsl_matrix *perm, gsl_matrix *gradR ) {
  computeGammaSr(R, myRorig, myTmpYr, true);

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
    computeGradFromYr(myTmpYr, myRorig, perm, gradR);
  }
}

void CostFunction::computeFuncAndPseudoJacobianLs( gsl_matrix *R, 
                       gsl_matrix *perm, gsl_vector *res, gsl_matrix *jac ) {
  computeGammaSr(R, myRorig, myTmpYr, true);
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
    computePseudoJacobianLsFromYr(myTmpYr, myRorig, perm, jac);
  } 
}

void CostFunction::computeCorrectionAndJacobian( gsl_matrix *R, 
                       gsl_matrix *perm, gsl_vector *res, gsl_matrix *jac ) {
  computeGammaSr(R, myRorig, myTmpYr, true);                  
  myGam->multInvGammaVector(myTmpYr);
  if (res != NULL) {
    gsl_vector_set_zero(res);
    myStruct->correctP(res, myRorig, myTmpYr, 1);
  }
  if (jac != NULL) {  
    computeJacobianOfCorrection(myTmpYr, R, perm, jac);
  } 
}

void CostFunction::computeCorrection( gsl_vector* p, gsl_matrix *R ) {
  try  {
    computeGammaSr(R, myRorig, myTmpYr, true);
  } catch (Exception *e) {
    if (!strncmp(e->getMessage(), "Gamma", 5)) {
      delete e;
      e = new Exception("Gamma matrix is singular. "
                        "Unable to compute the correction.\n");
    }
    throw e;
  }
  myGam->multInvGammaVector(myTmpYr);
  myStruct->correctP(p, myRorig, myTmpYr);
}

void CostFunction::computeDefaultRTheta( gsl_matrix *RTheta ) {
  int c_size1 = getN(), c_size2 = getNrow();
  size_t status = 0;
  size_t minus1 = -1;
  double tmp;

  gsl_matrix * tempc = gsl_matrix_alloc(c_size1, c_size2);
  if (myPhi == NULL) {
    gsl_matrix_memcpy(tempc, myMatr);
  } else {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myMatr, myPhi, 0, tempc);
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



