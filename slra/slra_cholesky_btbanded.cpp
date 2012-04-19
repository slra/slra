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

slraGammaCholeskyBTBanded::slraGammaCholeskyBTBanded( const slraStationaryStructure *s, 
    int r, double reg_gamma  ) :  slraGammaCholeskyBBanded(s, r, reg_gamma), myWs(s)  {
  myGamma = gsl_matrix_alloc(getD(), getD() * (getS() + 1));
  myWkTmp = gsl_matrix_alloc(getNplusD(), getD());
}  
  
slraGammaCholeskyBTBanded::~slraGammaCholeskyBTBanded() {
  gsl_matrix_free(myGamma);
  gsl_matrix_free(myWkTmp);
}

void slraGammaCholeskyBTBanded::computeGammak( gsl_matrix *R ) {
  size_t k;
  gsl_matrix_view submat;
  
  for (k = 0; k < getS(); k++) { /* compute brgamma_k = R' * w_k * R */
    submat = gsl_matrix_submatrix(myGamma, 0, k * getD(), getD(), getD());
    myWs->AtWkB(&submat.matrix, k, R, R, myWkTmp);
  }
  submat = gsl_matrix_submatrix(myGamma, 0, getS() * getD(), getD(), getD());
  gsl_matrix_set_zero(&submat.matrix);
}
  
void slraGammaCholeskyBTBanded::computeGammaUpperPart( gsl_matrix *R ) {
  computeGammak(R);
  
  int row_gam, col_gam, icor;
  double *gp = myPackedCholesky;
    
  for (int i = 0; i < d_times_s; i++) {
    for (int j = 0; j < getD(); j++) {
      icor = i + j + 1;
      gp[i + j * d_times_s] = gsl_matrix_get(myGamma, 
          icor % getD(), j + (getS() - (icor / getD())) * getD());
    }
  }
  for (int r = 1; r < getM(); r++) {
    gp +=  d_times_s * getD();
    memcpy(gp, myPackedCholesky, d_times_s * getD() * sizeof(double));
  }
}
  

/*void slraGammaCholeskyBBandedLH::computeCholeskyOfGamma( gsl_matrix *R ) { FAST version
  size_t i, j, sum_np, sum_nl, info = 0;
  gsl_vector_view R_row, TempR_row;
  
  for (i = 0; i < getM(); i++) {
    gsl_matrix blk_row = gsl_matrix_view_array_with_tda(myPackedCholesky, 
        (getS() + 1) * getD(), getD(), d_times_s_minus_1).matrix;
 
    gsl_matrix gamma_ij = gsl_matrix_submatrix(&blk_row, 0, 0, getD(), getD()).matrix;
    /* Compute Gamma_{i,i} * / 
    gsl_matrix_memcpy(myTempR, R);
    myWs->mulInvWij(myTempR, i);
    gsl_blas_dsyr2k(CblasLower, CblasTrans, 0.5, R, myTempR, 0.0, &gamma_ij);


    myWs->AtWijB(myTmpGammaij, j, i, R, R, myTempWktR);  
    
  
    /* Compute Gamma_{i, i+j} * / 
    for (j = 1; (j < getS()) && (j < getM() - i);  j++) {
      gamma_ij = gsl_matrix_submatrix(&blk_row, 0, j * getD(), getD(), getD()).matrix;
      myWs->AtWijB(&gamma_ij, i+j, i, R, R, myTempWktR);
    //  myWs->slraLayeredHankelStructure::AtWkB(&gamma_ij, -k, R, myTempR, myTempWktR);
    }
    if (getS() < getM() - i)  {
      gamma_ij = gsl_matrix_submatrix(&blk_row, 0, getS() * getD(), getD(), getD()).matrix;
      gsl_matrix_set_zero(&gamma_ij);
    }
  }

  dpbtrf_("U", &d_times_Mg, &d_times_s_minus_1, myPackedCholesky, &d_times_s, &info);


  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED * /
  }
}*/


void slraGammaCholeskySameDiagBTBanded::multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), myStruct->getMl(k) * myBase->getD());    
  
    myBase->multiplyInvPartCholeskyArray(sub_yr.vector.data, trans, 
        sub_yr.vector.size, sub_yr.vector.size);
  }
}

void slraGammaCholeskySameDiagBTBanded::multiplyInvGammaVector( gsl_vector * yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), myStruct->getMl(k) * myBase->getD());    
  
    myBase->multiplyInvPartGammaArray(sub_yr.vector.data, sub_yr.vector.size, sub_yr.vector.size);
  }
}
