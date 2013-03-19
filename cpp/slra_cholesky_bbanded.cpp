extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

SDependentCholesky::SDependentCholesky( const SDependentStructure *s,
     size_t D ) : myW(s), myD(D) {
  /* Calculate variables for FORTRAN routines */     
  s_minus_1 =  myW->getS() - 1;
  d_times_s =  myD * myW->getS();
  d_times_n = myW->getN() * myD;
  d_times_s_minus_1 = myD * myW->getS() - 1;
  
  /* Preallocate arrays */
  myPackedCholesky = (double*)malloc(d_times_n * d_times_s * sizeof(double));
  myTempR = gsl_matrix_alloc(myW->getM(), myD);
  myTempWktR = gsl_matrix_alloc(myW->getM(), myD);
  myTempGammaij = gsl_matrix_alloc(myD, myD);
}
  
SDependentCholesky::~SDependentCholesky() {
  free(myPackedCholesky);
  gsl_matrix_free(myTempR);
  gsl_matrix_free(myTempWktR);
  gsl_matrix_free(myTempGammaij);
}

void SDependentCholesky::multInvCholeskyVector( gsl_vector * yr, long trans ) {
  if (yr->stride != 1) {
    throw new Exception("Cannot mult vectors with stride != 1\n");
  }
  if (yr->size > d_times_n) {
    throw new Exception("yr->size > d_times_n\n");
  }
 
  size_t one = 1, info;
  dtbtrs_("U", (trans ? "T" : "N"), "N", &yr->size, &d_times_s_minus_1, &one, 
	        myPackedCholesky, &d_times_s, yr->data, &yr->size, &info);
}

void SDependentCholesky::multInvGammaVector( gsl_vector * yr ) {
  if (yr->stride != 1) {
    throw new Exception("Cannot mult vectors with stride != 1\n");
  }
  if (yr->size > d_times_n) {
    throw new Exception("yr->size > d_times_n\n");
  }

  size_t one = 1, info;
  dpbtrs_("U", &yr->size, &d_times_s_minus_1, &one, 
          myPackedCholesky, &d_times_s, yr->data, &yr->size, &info);  
}

void SDependentCholesky::calcGammaCholesky( const gsl_matrix *R, double reg_gamma ) {
  size_t info = 0;
  gsl_matrix m = gsl_matrix_view_array(myPackedCholesky, d_times_n, 
                     d_times_s).matrix;
  gsl_vector m_col = gsl_matrix_column(&m, d_times_s - 1).vector;
  
  computeGammaUpperPart(R);
  dpbtrf_("U", &d_times_n, &d_times_s_minus_1, myPackedCholesky, 
          &d_times_s, &info);
          
  if (info && reg_gamma > 0) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, 
        "Gamma matrix is singular (DPBTRF info = %d), "
        "adding regularization, reg = %f.\n", info, reg_gamma);
    computeGammaUpperPart(R, reg_gamma);

    dpbtrf_("U", &d_times_n, &d_times_s_minus_1,
        myPackedCholesky, &d_times_s, &info);
  }
  
  if (info) {
    throw new Exception("Gamma is singular (DPBTRF info = %d).\n", info); 
  }
}

void SDependentCholesky::computeGammaUpperPart( const gsl_matrix *R, double reg ) {
  gsl_matrix gamma_ij;
  gsl_vector diag;
  double *diagPtr =  myPackedCholesky;
  for (size_t i = 0; i < getN(); ++i, diagPtr += getS() * getD() * getD()) {
    if (getS() > 1) {
      gsl_matrix blk_row = 
          gsl_matrix_view_array_with_tda(diagPtr + d_times_s_minus_1, 
              (getS() + 1) * getD(), getD(), d_times_s_minus_1).matrix;
      for (size_t j = 0; (j <= getS()) && (j < getN() - i); j++) {
        gamma_ij = gsl_matrix_submatrix(&blk_row, j * getD(), 0, 
                                        getD(), getD()).matrix;
        if (j < getS()) {
          myW->AtWijB(myTempGammaij, i+j, i, R, R, myTempWktR);  
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
      myW->AtWijB(myTempGammaij, i, i, R, R, myTempWktR);  
      gamma_ij = gsl_matrix_view_array(diagPtr, getD(), getD()).matrix;
      shiftLowerTrg(&gamma_ij, myTempGammaij);
    }
  }
}
