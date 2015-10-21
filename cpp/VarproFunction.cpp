#include <memory.h>
#include <string.h>

#include "slra.h"

VarproFunction::VarproFunction( const gsl_vector *p, Structure *s, size_t d, 
                    gsl_matrix *Phi, bool isGCD ) : myP(NULL), myD(d), myStruct(s), 
                         myReggamma(SLRA_DEF_reggamma), myIsGCD(isGCD) {
  if (myStruct->getNp() > p->size) {
    throw new Exception("Inconsistent parameter vector\n");
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
    throw new Exception("Phi is not NULL\n");
  }

  if (d >= getNrow() || d <= 0) {
    throw new Exception("Incorrect rank given\n");
  }

  myP = gsl_vector_alloc(p->size);
  gsl_vector_memcpy(myP, p);
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
  myTmpJtJ = gsl_matrix_alloc(getNrow() * getD(), getNrow() * getD());
  myTmpEye = gsl_matrix_alloc(getNrow(), getNrow());
  gsl_matrix_set_identity(myTmpEye);
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
  gsl_vector_free(myPhiPermCol);
  gsl_matrix_free(myMatr);
  gsl_matrix_free(myTmpGradR);
  gsl_matrix_free(myTmpGradR2);
  gsl_vector_free(myTmpYr);
  gsl_matrix_free(myTmpJac);
  gsl_matrix_free(myTmpJac2);
  gsl_matrix_free(myTmpJtJ);
  gsl_matrix_free(myTmpEye);
  gsl_vector_free(myTmpJacobianCol);
  gsl_vector_free(myTmpCorr);
}

void VarproFunction::computeGammaSr( const gsl_matrix *Rt,
                                    gsl_vector *Sr, bool regularize_gamma ) {
  myGam->calcGammaCholesky(Rt, regularize_gamma ? myReggamma : 0);
  gsl_matrix SrMat = gsl_matrix_view_vector(Sr, getN(), getD()).matrix;
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myMatr, Rt, 0, &SrMat);
} 

void VarproFunction::fillZmatTmpJac( gsl_matrix *Zmatr, const gsl_vector* y,
                                     const gsl_matrix *Rt, double factor,
                                     int mult_gam ) {
  for (size_t j_1 = 0; j_1 < getM(); j_1++) {
    for (size_t i_1 = 0; i_1 < getD(); i_1++) {
      gsl_vector tJr = gsl_matrix_row(Zmatr, j_1 * getD() + i_1).vector;
    
      myDeriv->calcDijGammaYr(&tJr, Rt, j_1, i_1, y);
      gsl_vector_scale(&tJr, -factor);
      for (size_t k = 0; k < getN(); k++) {  /* Convert to vector strides */
        (*gsl_vector_ptr(&tJr, i_1 + k * getD())) +=
             gsl_matrix_get(myMatr, k, j_1);
      }
      
      if (mult_gam == 1) {
         myGam->multInvCholeskyVector(&tJr, 1);
      } else if (mult_gam == 2) {
         myGam->multInvGammaVector(&tJr);
      }
    }
  }
}



void VarproFunction::setPhiPermCol( size_t i, const gsl_matrix *perm,
                                    gsl_vector *phiPermCol ) {
  if (perm != NULL) {
    gsl_vector permCol = gsl_matrix_const_column(perm, i).vector;
    gsl_vector_memcpy(phiPermCol, &permCol);
  } else {
    gsl_vector phiCol = gsl_matrix_column(myTmpEye, i).vector;
    gsl_vector_memcpy(phiPermCol, &phiCol);
  }  
}

void VarproFunction::mulZmatPerm( gsl_vector* res, const gsl_matrix *Zmatr,
         const gsl_matrix *PsiT, size_t j_1, size_t i_1 ) {
  gsl_matrix subJ =
      gsl_matrix_view_array_with_tda(Zmatr->data + i_1 * Zmatr->tda, getM(),
                                     Zmatr->size2, getD() * Zmatr->tda).matrix;
  setPhiPermCol(j_1, PsiT, myPhiPermCol);
  gsl_blas_dgemv(CblasTrans, 1.0, &subJ, myPhiPermCol, 0.0, res); 
}

void VarproFunction::computePseudoJacobianLsFromYr( const gsl_vector* yr, 
         const gsl_matrix *Rt, const gsl_matrix *PsiT, gsl_matrix *pjac,
         double factor ) {
  fillZmatTmpJac(myTmpJac, yr, Rt, factor, 1);

  if (PsiT == NULL || PsiT->size1 == getNrow()) {
    size_t nrow = PsiT != NULL ? PsiT->size2 : getNrow();
    for (size_t j_1 = 0; j_1 < nrow; j_1++) {
      for (size_t i_1 = 0; i_1 < getD(); i_1++) {
        mulZmatPerm(myTmpJacobianCol, myTmpJac, PsiT, j_1, i_1);
        gsl_matrix_set_col(pjac, j_1 * getD() + i_1, myTmpJacobianCol);
      }
    }
  } else {
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, myTmpJac, PsiT, 0, pjac);
  }
}


void VarproFunction::computeJacobianOfCorrection( const gsl_vector* yr, 
         const gsl_matrix *Rt, const gsl_matrix *PsiT, gsl_matrix *jac ) {
  size_t nrow = PsiT != NULL ? PsiT->size2 : getNrow();
  gsl_matrix_set_zero(jac);
  gsl_matrix_set_zero(myTmpGradR);

  fillZmatTmpJac(myTmpJac, yr, Rt, 1, 2);

  if (PsiT == NULL || PsiT->size1 == getNrow()) {
    for (size_t j_1 = 0; j_1 < nrow; j_1++) {
      for (size_t i_1 = 0; i_1 < getD(); i_1++) {
        /* Compute first term (correction of Gam^{-1} z_{ij}) */
        gsl_vector_set_zero(myTmpCorr);
        mulZmatPerm(myTmpJacobianCol, myTmpJac, PsiT, j_1, i_1);
        myStruct->multByGtUnweighted(myTmpCorr, Rt, myTmpJacobianCol, -1, 1);

        /* Compute second term (gamma * dG_{ij} * yr) */
        gsl_matrix_set_zero(myTmpGradR);
        setPhiPermCol(j_1, PsiT, myPhiPermCol);
        gsl_matrix_set_col(myTmpGradR, i_1, myPhiPermCol);
        myStruct->multByGtUnweighted(myTmpCorr, myTmpGradR, yr, -1, 1);

        myStruct->multByWInv(myTmpCorr, 1);

        gsl_vector jac_col = gsl_matrix_column(jac, j_1 * getD() + i_1).vector;
        gsl_vector_memcpy(&jac_col, myTmpCorr);
      }
    }
  } else {
    for (size_t i = 0; i < PsiT->size2; i++) {
      gsl_vector PsiRow = gsl_matrix_const_column(PsiT, i).vector;
      
      /* Compute first term (correction of Gam^{-1} z_{ij}) */
      gsl_vector_set_zero(myTmpCorr);
      gsl_blas_dgemv(CblasTrans, 1.0, myTmpJac, &PsiRow, 0.0, myTmpJacobianCol);
      myStruct->multByGtUnweighted(myTmpCorr, Rt, myTmpJacobianCol, -1, 1);

      /* Compute second term (gamma * dG_{ij} * yr) */
      //TODO: relies on the fact that tda = size2
      gsl_vector GradVec = gsl_vector_view_array(myTmpGradR->data,
                                myTmpGradR->size1 * myTmpGradR->tda).vector;
      gsl_vector_memcpy(&GradVec, &PsiRow);
      myStruct->multByGtUnweighted(myTmpCorr, myTmpGradR, yr, -1, 1);
      
      myStruct->multByWInv(myTmpCorr, 1);
      
      gsl_vector jac_col = gsl_matrix_column(jac, i).vector;
      gsl_vector_memcpy(&jac_col, myTmpCorr);
    }
  }
}

void VarproFunction::computeGradFromYr( const gsl_vector* yr, 
         const gsl_matrix *Rt, const gsl_matrix *perm, gsl_matrix *gradR ) {
  gsl_matrix_const_view yr_matr = gsl_matrix_const_view_vector(yr, getN(), getD());
  myDeriv->calcYtDgammaY(myTmpGradR, Rt, &yr_matr.matrix);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myMatr, &yr_matr.matrix,
                -1.0, myTmpGradR);
  gsl_matrix_memcpy(myTmpGradR2, myTmpGradR);

  if (perm != NULL) {
    if (perm->size1 == getNrow() && perm->size2 == gradR->size1) { // TODO: improve on this
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, perm, myTmpGradR2, 0.0, gradR);
    } else {
      gsl_vector vecTmpGradR2 = gsl_vector_view_array(myTmpGradR2->data,
                                     myTmpGradR2->size1 * myTmpGradR2->size2).vector, 
                 vecGradR =  gsl_vector_view_array(gradR->data,
                                     gradR->size1 * gradR->size2).vector;
      gsl_blas_dgemv(CblasTrans, 1.0, perm, &vecTmpGradR2, 0.0, &vecGradR); 
    }
  } else {
    gsl_matrix_memcpy(gradR, myTmpGradR2);
  }
}

void VarproFunction::computeFuncAndGrad( const gsl_matrix* Rt, double * f,
                                       const gsl_matrix *perm, gsl_matrix *gradR ) {
  computeGammaSr(Rt, myTmpYr, true);

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
    computeGradFromYr(myTmpYr, Rt, perm, gradR);
    if (myIsGCD) {
      gsl_matrix_scale(gradR, -1);
    }
  }
}

void VarproFunction::computeCorrectionAndJacobian( const gsl_matrix *Rt,
         const gsl_matrix *perm, gsl_vector *res, gsl_matrix *jac ) {
  computeGammaSr(Rt, myTmpYr, true);
  myGam->multInvGammaVector(myTmpYr);
  if (res != NULL) {
    if (myIsGCD) {
      gsl_vector_memcpy(res, getP());
      myStruct->multByGtUnweighted(res, Rt, myTmpYr, -1, 1);
    } else {
      gsl_vector_set_zero(res);
      myStruct->multByGtUnweighted(res, Rt, myTmpYr, -1, 1);
    }
    myStruct->multByWInv(res, 1);
  }
  if (jac != NULL) {  
    computeJacobianOfCorrection(myTmpYr, Rt, perm, jac);
  } 
}

void VarproFunction::computePhat( gsl_vector* p, const gsl_matrix *Rt ) {
  try  {
    computeGammaSr(Rt, myTmpYr, true);
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
    myStruct->multByGtUnweighted(p, Rt, myTmpYr, 1, 1, false);
  } else {
    myStruct->multByGtUnweighted(p, Rt, myTmpYr, -1, 1);
    myStruct->multByWInv(p, 2);
    gsl_vector_add(p, getP());
  }
}

void VarproFunction::computeFuncAndPseudoJacobianLs( const gsl_matrix *Rt,
         gsl_matrix *perm, gsl_vector *res, gsl_matrix *jac, double factor ) {
  if (myIsGCD)  {
    throw new Exception("Pseudojacobian not allowed for GCD computations\n");
  }
  computeGammaSr(Rt, myTmpYr, true);
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
    computePseudoJacobianLsFromYr(myTmpYr, Rt, perm, jac, factor);
  } 
}

void VarproFunction::computeDefaultRTheta( gsl_matrix *RTheta ) {
  size_t c_size1 = getN(), c_size2 = getNrow();
  size_t status = 0;
  size_t minus1 = -1;
  double tmp;

  gsl_matrix * tempc = gsl_matrix_alloc(c_size1, c_size2);
  gsl_matrix_memcpy(tempc, myMatr);
  
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

void VarproFunction::computeJtJmulE( const gsl_matrix* R, const gsl_matrix* E, gsl_matrix *out, int useJtJ ) {
  gsl_vector vecE = gsl_vector_const_view_array(E->data, E->size1 * E->size2).vector;   
  gsl_vector vecOut = gsl_vector_const_view_array(out->data, out->size1 * out->size2).vector;   
  
  if (R->size1 * R->size2 != 0) { 
    computeFuncAndPseudoJacobianLs(R, myTmpEye, NULL, myTmpJac2);
  } /* Otherwise use the precomputed Jacobian */
  if (useJtJ) {
    if (R->size1 * R->size2 != 0) { 
	  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, myTmpJac2, myTmpJac2, 0, myTmpJtJ);
	} /* Otherwise use the precomputed JtJ */
    gsl_blas_dgemv(CblasNoTrans, 2.0, myTmpJtJ, &vecE, 0.0, &vecOut);  
  } else {
    gsl_blas_dgemv(CblasNoTrans, 1.0, myTmpJac2, &vecE, 0.0, myTmpJacobianCol);  
    gsl_blas_dgemv(CblasTrans, 2.0, myTmpJac2, myTmpJacobianCol, 0.0, &vecOut);  
  }
}



