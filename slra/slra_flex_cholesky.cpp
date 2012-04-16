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


slraBTBGammaCholesky::slraBTBGammaCholesky( const slraWkInterface *s, int r, int Mg, 
     int use_slicot, double reg_gamma  ) : 
     myMg(Mg), myN(r), myD(s->getNplusD()-r), my_use_slicot(use_slicot), my_reg_gamma(reg_gamma)  {
     
  myW = s;
    
  /* Preallocate arrays */
  myGammaVec = (double*) malloc(myD * myD * (myW->getS() + 1) * sizeof(double));
  myGamma = gsl_matrix_alloc(myD, myD * (myW->getS() + 1));
  myWkTmp = gsl_matrix_alloc(myD, myN + myD);
  myPackedCholesky = (double *)malloc(myMg * myD * myD * myW->getS() * 
                              sizeof(double));
  myCholeskyWorkSize = 1 + myW->getS() * myD * myD + /* pDW */ 
                       3 * myD + /* 3 * K */
                       mymax(myW->getS(), myMg - myW->getS()) * myD * myD;
  myCholeskyWork = (double *)malloc(myCholeskyWorkSize * sizeof(double));                       
                       

                       
  /* Calculate variables for FORTRAN routines */     
  s_minus_1 = myW->getS() - 1;
  d_times_s = myD * myW->getS();
  d_times_Mg = myMg * myD;
  d_times_s_minus_1 = myD * myW->getS() - 1;
}
  
slraBTBGammaCholesky::~slraBTBGammaCholesky() {
  free(myGammaVec);
  gsl_matrix_free(myGamma);
  gsl_matrix_free(myWkTmp);
  free(myPackedCholesky);
  free(myCholeskyWork);
}
  
void slraBTBGammaCholesky::computeCholeskyOfGamma( gsl_matrix *R )  {
  int k;
  size_t info = 0;
  gsl_matrix_view submat;
  const size_t zero = 0;

  /* compute brgamma_k = R' * w_k * R */
  for (k = 0; k < myW->getS(); k++) {
    submat = gsl_matrix_submatrix(myGamma, 0, k * myD, myD, myD);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, R, myW->getWk(k), 0.0, myWkTmp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myWkTmp, R, 0.0, 
                   &submat.matrix);
  }
  submat = gsl_matrix_submatrix(myGamma, 0, myW->getS() * myD, myD, myD);
  gsl_matrix_set_zero(&submat.matrix);
    
#ifdef USE_SLICOT    
  if (my_use_slicot) { /* use SLICOT */
    gsl_matrix_vectorize(myGammaVec, myGamma);
    /* Cholesky factorization of Gamma */
    mb02gd_("R", "N", &myD, &myMg, &s_minus_1, &zero, 
            &myMg, myGammaVec, &myD, myPackedCholesky, &d_times_s, 
            myCholeskyWork, &myCholeskyWorkSize, &info); /**/
  } else { /* use LAPACK */
#endif  
    int i, j, r, row_gam, col_gam, icor;
    double *gp = myPackedCholesky;
    
    for (i = 0; i < d_times_s; i++) {
      for (j = 0; j < myD; j++) {
        icor = i + j + 1;
        gp[i + j * d_times_s] = gsl_matrix_get(myGamma, 
            icor % myD, j + (myW->getS() - (icor / myD)) * myD);
      }
    }
    
    for (r = 1; r < myMg; r++) {
      gp +=  d_times_s * myD;
      memcpy(gp, myPackedCholesky, d_times_s * myD * sizeof(double));
    }
    
    dpbtrf_("U", &d_times_Mg, &d_times_s_minus_1, myPackedCholesky, 
            &d_times_s, &info);
#ifdef USE_SLICOT
  }
#endif  
  
  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED */
  }
}
  
void slraBTBGammaCholesky::multiplyInvPartCholeskyArray( double * yr, int trans, size_t size, size_t chol_size ) {
  size_t info;
  size_t total_cols = size / chol_size;

  dtbtrs_("U", (trans ? "T" : "N"), "N", 
          &chol_size, &d_times_s_minus_1, &total_cols, 
	  myPackedCholesky, &d_times_s, yr, &chol_size, &info);
}
  
void slraBTBGammaCholesky::multiplyInvPartGammaArray( double * yr, size_t size, size_t chol_size ) {
  size_t info;
  size_t total_cols = size / chol_size; 
  
  dpbtrs_("U", &chol_size, &d_times_s_minus_1, &total_cols, 
          myPackedCholesky, &d_times_s, yr, &chol_size, &info);  
}

void slraSameDiagBTBGammaCholesky::multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), myStruct->getMl(k) * myBase->getD());    
  
    myBase->multiplyInvPartCholeskyArray(sub_yr.vector.data, trans, 
        sub_yr.vector.size, sub_yr.vector.size);
  }
}

void slraSameDiagBTBGammaCholesky::multiplyInvGammaVector( gsl_vector * yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), myStruct->getMl(k) * myBase->getD());    
  
    myBase->multiplyInvPartGammaArray(sub_yr.vector.data, sub_yr.vector.size, sub_yr.vector.size);
  }
}

typedef slraGammaCholesky* pslraGammaCholesky;

slraDiagGammaCholesky::slraDiagGammaCholesky( const slraStripedStructure *s, int r, double reg_gamma  ) :
    myStruct(s) {
  myD = s->getNplusD() - r;
  myGamma = new pslraGammaCholesky[myStruct->getBlocksN()];
  for (int k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k] = myStruct->getBlock(k)->createGammaComputations(r, reg_gamma);
  }
}    

slraDiagGammaCholesky::~slraDiagGammaCholesky() {
  if (myGamma != NULL) {
    for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
      if (myGamma[k] != NULL) {
        delete myGamma[k];
      }
    }
    delete[] myGamma;
  }
}

void slraDiagGammaCholesky::multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
  
    myGamma[k]->multiplyInvCholeskyVector(&sub_yr.vector, trans);
  }
}

void slraDiagGammaCholesky::multiplyInvGammaVector( gsl_vector * yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
  
    myGamma[k]->multiplyInvGammaVector(&sub_yr.vector);
  }
}


void slraDiagGammaCholesky::computeCholeskyOfGamma( gsl_matrix *R ) {
  for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k]->computeCholeskyOfGamma(R);  
  }
}




