#include "slra.h"

MuDependentCholesky::MuDependentCholesky( const MuDependentStructure *s,
                                        size_t D ) : myW(s), myD(D) {
  /* Calculate variables for FORTRAN routines */     
  myMu_1 =  myW->getMu() - 1;
  myDMu =  myD * myW->getMu();
  myDN = myW->getN() * myD;
  myDMu_1 = myD * myW->getMu() - 1;
  /* Preallocate arrays */
  myPackedCholesky = (double*)malloc(myDN * myDMu * sizeof(double));
  myTempR = gsl_matrix_alloc(myW->getM(), myD);
  myTempVktR = gsl_matrix_alloc(myW->getM(), myD);
  myTempGammaij = gsl_matrix_alloc(myD, myD);
}
  
MuDependentCholesky::~MuDependentCholesky() {
  free(myPackedCholesky);
  gsl_matrix_free(myTempR);
  gsl_matrix_free(myTempVktR);
  gsl_matrix_free(myTempGammaij);
}

void MuDependentCholesky::multInvCholeskyVector( gsl_vector * yr, long trans ) {
  if (yr->stride != 1) {
    throw new Exception("Cannot mult vectors with stride != 1\n");
  }
  if (yr->size > myDN) {
    throw new Exception("yr->size > myDN\n");
  }
  size_t one = 1, info;
  dtbtrs_("U", (trans ? "T" : "N"), "N", &yr->size, &myDMu_1, &one, 
	        myPackedCholesky, &myDMu, yr->data, &yr->size, &info);
}

void MuDependentCholesky::multInvGammaVector( gsl_vector * yr ) {
  if (yr->stride != 1) {
    throw new Exception("Cannot mult vectors with stride != 1\n");
  }
  if (yr->size > myDN) {
    throw new Exception("yr->size > myDN\n");
  }
  size_t one = 1, info;
  dpbtrs_("U", &yr->size, &myDMu_1, &one, 
          myPackedCholesky, &myDMu, yr->data, &yr->size, &info);  
}

void MuDependentCholesky::calcGammaCholesky( const gsl_matrix *R, double reg ) {
  size_t info = 0;
  gsl_matrix m = gsl_matrix_view_array(myPackedCholesky, myDN, myDMu).matrix;
  gsl_vector m_col = gsl_matrix_column(&m, myDMu - 1).vector;
  
  computeGammaUpperPart(R);
  dpbtrf_("U", &myDN, &myDMu_1, myPackedCholesky, &myDMu, &info);
          
  if (info && reg > 0) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, "Gamma is singular (DPBTRF info = %d), "
        "adding regularization, reg = %f.\n", info, reg);
    computeGammaUpperPart(R, reg);
    dpbtrf_("U", &myDN, &myDMu_1, myPackedCholesky, &myDMu, &info);
  }
  
  if (info) {
    throw new Exception("Gamma is singular (DPBTRF info = %d).\n", info); 
  }
}

void MuDependentCholesky::computeGammaUpperPart( const gsl_matrix *R, double reg ) {
  gsl_matrix gamma_ij;
  gsl_vector diag;
  double *diagPtr =  myPackedCholesky;
  for (size_t i = 0; i < getN(); ++i, diagPtr += getMu() * getD() * getD()) {
    if (getMu() > 1) {
      gsl_matrix blk_row = 
          gsl_matrix_view_array_with_tda(diagPtr + myDMu_1, 
              (getMu() + 1) * getD(), getD(), myDMu_1).matrix;
      for (size_t j = 0; (j <= getMu()) && (j < getN() - i); j++) {
        gamma_ij = gsl_matrix_submatrix(&blk_row, j * getD(), 0, 
                                        getD(), getD()).matrix;
        if (j < getMu()) {
          myW->AtVijB(myTempGammaij, i+j, i, R, R, myTempVktR);  
        } else {
          gsl_matrix_set_zero(myTempGammaij);
        }
        if (j == 0) {
          if (reg > 0) {
            diag = gsl_matrix_diagonal(myTempGammaij).vector;
            gsl_vector_add_constant(&diag, reg);
          }
        
          copyLowerTrg(&gamma_ij, myTempGammaij);
        } else {
          gsl_matrix_memcpy(&gamma_ij, myTempGammaij);
        }
      }
    } else {
      myW->AtVijB(myTempGammaij, i, i, R, R, myTempVktR);  
      gamma_ij = gsl_matrix_view_array(diagPtr, getD(), getD()).matrix;
      shiftLowerTrg(&gamma_ij, myTempGammaij);
    }
  }
}
