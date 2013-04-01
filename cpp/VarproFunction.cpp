#include <memory.h>
#include <string.h>

#include "slra.h"

VarproFunction::VarproFunction( const gsl_vector *p, Structure *s, size_t d, 
                    gsl_matrix *Phi, bool isGCD ) : myP(NULL), myD(d), myStruct(s), 
                         myReggamma(SLRA_DEF_reggamma), myIsGCD(isGCD) {
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

  if (myStruct->getN() * getD() * myStruct->getM() * getD() >= 10000000L) {
    throw new Exception("Too much memory required: the Jacobian would have "
                  "more than 10^8 elements. This is currently not allowed.\n");
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
  myTmpJac2 = gsl_matrix_alloc(myStruct->getN() * getD(), getNrow() * getD());
  myTmpCorr = gsl_vector_alloc(myStruct->getNp());
  if (myIsGCD) {
    gsl_vector_memcpy(myTmpCorr, getP());
    myStruct->multByWInv(myTmpCorr, 1);
    myStruct->fillMatrixFromP(myMatr, myTmpCorr);
    gsl_blas_ddot(myTmpCorr, myTmpCorr, &myPWnorm2);
  } else {
    myStruct->fillMatrixFromP(myMatr, getP());
  }
}
  
VarproFunction::~VarproFunction() {
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
  gsl_matrix_free(myTmpJac2);
  gsl_vector_free(myTmpJacobianCol);
  gsl_vector_free(myTmpCorr);
}

void VarproFunction::computeGammaSr( const gsl_matrix *R, 
                      gsl_matrix *Rorig, gsl_vector *Sr, bool regularize_gamma ) {
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPhi, R, 0, Rorig);
  myGam->calcGammaCholesky(Rorig, regularize_gamma ? myReggamma : 0);
  gsl_matrix_view SrMat = gsl_matrix_view_vector(Sr, getN(), getD()); 
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myMatr, Rorig, 0, &SrMat.matrix);
} 

void VarproFunction::fillZmatTmpJac( const gsl_vector* yr, 
                                        const gsl_matrix *Rorig, double factor ) {
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

void VarproFunction::setPhiPermCol( size_t i, const gsl_matrix *perm ) {
  if (perm != NULL) {
    gsl_vector permCol = gsl_matrix_const_column(perm, i).vector;
    gsl_blas_dgemv(CblasNoTrans, 1.0, myPhi, &permCol, 0.0, myPhiPermCol);
  } else {
    gsl_vector phiCol = gsl_matrix_column(myPhi, i).vector;
    gsl_vector_memcpy(myPhiPermCol, &phiCol);
  }  
}

void VarproFunction::mulZmatPerm( gsl_vector* res, const gsl_matrix *perm,
                                  size_t i, size_t j ) {
  gsl_matrix subJ = gsl_matrix_submatrix(myTmpJac, j * getM(), 0, getM(),
                         myTmpJac->size2).matrix;
  setPhiPermCol(i, perm);    
  gsl_blas_dgemv(CblasTrans, 1.0, &subJ, myPhiPermCol, 0.0, res); 
}

void VarproFunction::computePseudoJacobianLsFromYr( const gsl_vector* yr, 
         const gsl_matrix *Rorig, const gsl_matrix *perm, gsl_matrix *pjac, 
         double factor ) {
  size_t nrow = perm != NULL ? perm->size2 : getNrow();
  fillZmatTmpJac(yr, Rorig, factor);
  
  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < getD(); j++) {
      mulZmatPerm(myTmpJacobianCol, perm, i, j);
      myGam->multInvCholeskyVector(myTmpJacobianCol, 1);
      gsl_matrix_set_col(pjac, i * getD() + j, myTmpJacobianCol);  
    }
  }
}


void VarproFunction::computeJacobianOfCorrection( const gsl_vector* yr, 
         const gsl_matrix *Rorig, const gsl_matrix *perm, gsl_matrix *jac ) {
  size_t nrow = perm != NULL ? perm->size2 : getNrow();
  gsl_matrix_set_zero(jac);
  gsl_matrix_set_zero(myTmpGradR);

  fillZmatTmpJac(yr, Rorig, 1);

  for (size_t i = 0; i < perm->size2; i++) {
    for (size_t j = 0; j < getD(); j++) {
      gsl_vector_view jac_col = gsl_matrix_column(jac, i * getD() + j);

      /* Compute first term (correction of Gam^{-1} z_{ij}) */
      gsl_vector_set_zero(myTmpCorr);
      mulZmatPerm(myTmpJacobianCol, perm, i, j);
      myGam->multInvGammaVector(myTmpJacobianCol);
      myStruct->multByGtUnweighted(myTmpCorr, Rorig, myTmpJacobianCol, -1, 1);

      /* Compute second term (gamma * dG_{ij} * yr) */
      gsl_matrix_set_zero(myTmpGradR);
      setPhiPermCol(i, perm);
      gsl_matrix_set_col(myTmpGradR, j, myPhiPermCol);
      myStruct->multByGtUnweighted(myTmpCorr, myTmpGradR, yr, -1, 1);

      myStruct->multByWInv(myTmpCorr, 1);
      gsl_vector_memcpy(&jac_col.vector, myTmpCorr);
    }
  }
}

void VarproFunction::computeGradFromYr( const gsl_vector* yr, 
         const gsl_matrix *Rorig, const gsl_matrix *perm, gsl_matrix *gradR ) {
  gsl_matrix_const_view yr_matr = gsl_matrix_const_view_vector(yr, getN(), getD());
  myDeriv->calcYrtDgammaYr(myTmpGradR, Rorig, yr);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myMatr, &yr_matr.matrix,
                -1.0, myTmpGradR);
  if (myPhi != NULL) {
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, myPhi, myTmpGradR,
                 0.0, myTmpGradR2);
  } else {
    gsl_matrix_memcpy(myTmpGradR2, myTmpGradR);
  }
  if (perm != NULL) {
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, perm, myTmpGradR2, 0.0, gradR);
  } else {
    gsl_matrix_memcpy(gradR, myTmpGradR2);
  }
}

void VarproFunction::computeFuncAndGrad( const gsl_matrix* R, double * f, 
                                       const gsl_matrix *perm, gsl_matrix *gradR ) {
  computeGammaSr(R, myRorig, myTmpYr, true);

  if (f != NULL) {
    myGam->multInvCholeskyVector(myTmpYr, 1);
    gsl_blas_ddot(myTmpYr, myTmpYr, f);
    if (myIsGCD) {
      *f =  myPWnorm2 - *f;
    }
  }
  if (gradR != NULL) {
    if (f != NULL) {
      myGam->multInvCholeskyVector(myTmpYr, 0);
    } else {
      myGam->multInvGammaVector(myTmpYr);
    }
    computeGradFromYr(myTmpYr, myRorig, perm, gradR);
    if (myIsGCD) {
      gsl_matrix_scale(gradR, -1);
    }
  }
}

void VarproFunction::computeCorrectionAndJacobian( const gsl_matrix *R, 
         const gsl_matrix *perm, gsl_vector *res, gsl_matrix *jac ) {
  computeGammaSr(R, myRorig, myTmpYr, true);                  
  myGam->multInvGammaVector(myTmpYr);
  if (res != NULL) {
    if (myIsGCD) {
      gsl_vector_memcpy(res, getP());
      myStruct->multByGtUnweighted(res, myRorig, myTmpYr, -1, 1);
    } else {
      gsl_vector_set_zero(res);
      myStruct->multByGtUnweighted(res, myRorig, myTmpYr, -1, 1);
    }
    myStruct->multByWInv(res, 1);
  }
  if (jac != NULL) {  
    computeJacobianOfCorrection(myTmpYr, R, perm, jac);
  } 
}

void VarproFunction::computePhat( gsl_vector* p, const gsl_matrix *R ) {
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
  
  gsl_vector_set_zero(p);
  if (myIsGCD) {
    myStruct->multByGtUnweighted(p, myRorig, myTmpYr, 1, 1, false);
  } else {
    myStruct->multByGtUnweighted(p, myRorig, myTmpYr, -1, 1);
    myStruct->multByWInv(p, 2);
    gsl_vector_add(p, getP());
  }
}


void VarproFunction::computeFuncAndPseudoJacobianLs( const gsl_matrix *R, 
         gsl_matrix *perm, gsl_vector *res, gsl_matrix *jac, double factor ) {
  if (myIsGCD)  {
    throw new Exception("Pseudojacobian not allowed for GCD computations\n");
  }
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
    computePseudoJacobianLsFromYr(myTmpYr, myRorig, perm, jac, factor);
  } 
}

void VarproFunction::computeDefaultRTheta( gsl_matrix *RTheta ) {
  size_t c_size1 = getN(), c_size2 = getNrow();
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



