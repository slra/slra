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

slraGammaCholeskyBTBanded::slraGammaCholeskyBTBanded( const slraStationaryStructure *s, int r, double reg_gamma  ) : 
    slraGammaCholeskyBBanded(s, r, reg_gamma), myWs(s)  {
  myGamma = gsl_matrix_alloc(getD(), getD() * (getS() + 1));
  myWkTmp = gsl_matrix_alloc(getNplusD(), getD());
#ifdef USE_SLICOT
  myGammaVec = (double*) malloc(getD() * getD() * (getS() + 1) * sizeof(double));
  myCholeskyWorkSize = 1 + getS() * getD() * getD() + /* pDW */ 
                       3 * getD() + /* 3 * K */
                       mymax(getS(), getM() - getS()) * getD() * getD();
  myCholeskyWork = (double *)malloc(myCholeskyWorkSize * sizeof(double));                       
#endif  
}
  
slraGammaCholeskyBTBanded::~slraGammaCholeskyBTBanded() {
  gsl_matrix_free(myGamma);
  gsl_matrix_free(myWkTmp);
#ifdef USE_SLICOT  
  free(myGammaVec);
  free(myCholeskyWork);
#endif  
}
  
void slraGammaCholeskyBTBanded::computeCholeskyOfGamma( gsl_matrix *R )  {
  size_t info = 0, k;
  gsl_matrix_view submat;
  const size_t zero = 0;

  for (k = 0; k < getS(); k++) { /* compute brgamma_k = R' * w_k * R */
    submat = gsl_matrix_submatrix(myGamma, 0, k * getD(), getD(), getD());
    myWs->AtWkB(&submat.matrix, k, R, R, myWkTmp);
  }
  submat = gsl_matrix_submatrix(myGamma, 0, getS() * getD(), getD(), getD());
  gsl_matrix_set_zero(&submat.matrix);
    
#ifdef USE_SLICOT    
  size_t D = getD(), Mg = getM();
  gsl_matrix_vectorize(myGammaVec, myGamma);
    
  mb02gd_("R", "N", &D, &Mg, &s_minus_1, &zero, 
          &myMg, myGammaVec, &D, myPackedCholesky, &d_times_s, 
          myCholeskyWork, &myCholeskyWorkSize, &info); /**/
#else /* USE_SLICOT */
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
  
  dpbtrf_("U", &d_times_Mg, &d_times_s_minus_1, myPackedCholesky, 
          &d_times_s, &info);
#endif /* USE_SLICOT */       

  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED */
  }
}


void slraGammaCholeskyBBandedLH::computeCholeskyOfGamma( gsl_matrix *R ) {
  size_t i, l, k, j, sum_np, sum_nl, info = 0;
  gsl_vector_view R_row, TempR_row;
  
  for (i = 0; i < getM(); i++) {
    for (sum_np = i, sum_nl = 0, l = 0; l < myWs->getQ(); 
         sum_np += myWs->getLayerNp(k), sum_nl += myWs->getLayerLag(l), ++l) {
      for (k = 0; k < myWs->getLayerLag(l); k++) {
        double w = myWs->getInvSqrtWeights(sum_np + k);
        TempR_row = gsl_matrix_row(myTempR, sum_nl + k);
        gsl_vector_set_all(&TempR_row.vector, w * w);

        R_row = gsl_matrix_row(R, sum_nl + k);
        gsl_vector_mul(&TempR_row.vector, &R_row.vector);
      }
    }

    gsl_matrix_view blk_row = gsl_matrix_view_array_with_tda(myPackedCholesky, 
        (getS() + 1) * getD(), getD(), d_times_s_minus_1);
 
    /* Compute Gamma_{i,i} */ 
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, R, myTempR, 0.0, myTempGammaij);
    gsl_matrix_view gamma_ij = gsl_matrix_submatrix(&blk_row.matrix, 0, 0, getD(), getD());
    for (k = 0; k < getD(); k++) {
      for (l = 0; l < k; l++) {
        gsl_matrix_set(&gamma_ij.matrix, l, k, gsl_matrix_get(myTempGammaij, l, k));
      }
    }
  
    /* Compute Gamma_{i, i+j} */ 
    for (j = 1; j < mymin(getS() + 1, getM() - i);  j++) {
      gamma_ij = gsl_matrix_submatrix(&blk_row.matrix, 0, j * getD(), getD(), getD());
      if (j < getS()) {
        myWs->AtWkB(&gamma_ij.matrix, -k, R, myTempR, myTempWktR);
      } else {
        gsl_matrix_set_zero(myTempGammaij);
      }
    }
  }

  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED */
  }
}



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
