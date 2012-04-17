#include <memory.h>
#include <cstdarg>

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

slraGammaCholeskyBBanded::slraGammaCholeskyBBanded( const slraStructure *s, int r, double reg_gamma ) :
     myW(s), myD(s->getNplusD()-r), my_reg_gamma(reg_gamma) {
  /* Preallocate arrays */
  myPackedCholesky = (double *)malloc(myW->getM() * myD * myD * myW->getS() * sizeof(double));

  /* Calculate variables for FORTRAN routines */     
  s_minus_1 =  myW->getS() - 1;
  d_times_s =  myD * myW->getS();
  d_times_Mg = myW->getM() * myD;
  d_times_s_minus_1 = myD * myW->getS() - 1;
}
  
slraGammaCholeskyBBanded::~slraGammaCholeskyBBanded() {
  free(myPackedCholesky);
}

void slraGammaCholeskyBBanded::multiplyInvPartCholeskyArray( double * yr, int trans, size_t size, size_t chol_size ) {
  size_t info, total_cols = size / chol_size;

  dtbtrs_("U", (trans ? "T" : "N"), "N", 
          &chol_size, &d_times_s_minus_1, &total_cols, 
	  myPackedCholesky, &d_times_s, yr, &chol_size, &info);
}
  
void slraGammaCholeskyBBanded::multiplyInvPartGammaArray( double * yr, size_t size, size_t chol_size ) {
  size_t info, total_cols = size / chol_size; 
  
  dpbtrs_("U", &chol_size, &d_times_s_minus_1, &total_cols, 
          myPackedCholesky, &d_times_s, yr, &chol_size, &info);  
}

void slraGammaCholeskyBBanded::multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
  if (yr->stride != 1) {
    throw new slraException("Cannot multiply vectors with stride != 1\n");
  }
  multiplyInvPartCholeskyArray(yr->data, trans, yr->size, d_times_Mg);
}

void slraGammaCholeskyBBanded::multiplyInvGammaVector( gsl_vector * yr ) {
  if (yr->stride != 1) {
    throw new slraException("Cannot multiply vectors with stride != 1\n");
  }
  multiplyInvPartGammaArray(yr->data, yr->size, d_times_Mg);
}

void slraGammaCholeskyBBanded::multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans ) {
  if (yr_matr->size2 != yr_matr->tda) {
    slraGammaCholesky::multiplyInvCholeskyTransMatrix(yr_matr, trans);
  } else {
    multiplyInvPartCholeskyArray(yr_matr->data, trans, yr_matr->size1 * yr_matr->size2, d_times_Mg);
  }
}

