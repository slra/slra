#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

slraDGammaBTBanded::slraDGammaBTBanded( const slraStationaryStructure *s, int r ) :
    myD(s->getNplusD() - r), myW(s) {
  
  myTempWkColRow = gsl_vector_alloc(myW->getNplusD());
  myDGammaVec = gsl_vector_alloc(myD * (2 * myW->getS() - 1));
  myDGammaTrMat = gsl_matrix_alloc(myD, 2 * myW->getS() - 1);
  myDGamma = gsl_matrix_alloc(myD, myD * (2 * myW->getS() - 1));
  myTmpCol = gsl_vector_alloc(myW->getM());
  myWk_R =  gsl_matrix_alloc(myW->getNplusD(), myD);
  myWkT_R = gsl_matrix_alloc(myW->getNplusD(), myD);
  myN_k = gsl_matrix_alloc(myD, myD);
}

slraDGammaBTBanded::~slraDGammaBTBanded() {
  gsl_vector_free(myTempWkColRow);
  gsl_vector_free(myTmpCol);
  gsl_vector_free(myDGammaVec);
  gsl_matrix_free(myDGammaTrMat);
  gsl_matrix_free(myDGamma);
  gsl_matrix_free(myWk_R);
  gsl_matrix_free(myWkT_R);
  gsl_matrix_free(myN_k);
}

void slraDGammaBTBanded::calcYrtDgammaYr( gsl_matrix *mgrad_r, 
         gsl_matrix *R, gsl_vector *yr ) {
  int m = yr->size / myD;
  gsl_matrix Yr = gsl_matrix_view_vector(yr, m, myD).matrix, YrL, YrR;

  gsl_matrix_set_zero(mgrad_r);
  for (int k = 0; k < myW->getS(); k++) {
    YrL = gsl_matrix_submatrix(&Yr, 0, 0, m - k, myD).matrix;
    YrR = gsl_matrix_submatrix(&Yr, k, 0, m - k, myD).matrix;
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &YrL, &YrR, 0.0, myN_k);

    myW->WkB(myWk_R, k, R);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 2.0, myWk_R, myN_k, 1.0, mgrad_r);

    if (k > 0) {
      myW->WkB(myWk_R, -k, R);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 2.0, myWkT_R, myN_k, 1.0, mgrad_r);
    }
  }       
}

void slraDGammaBTBanded::calcDijGammaYr( gsl_vector *res, 
         gsl_matrix *R, gsl_matrix *perm, int i, int j,  gsl_vector *yr ) {
  gsl_vector gv_sub, perm_col = gsl_matrix_column(perm, i).vector, dgamma_j_row;

  for (int k = 1 - myW->getS(); k < myW->getS(); k++) {
    gv_sub = gsl_vector_subvector(myDGammaVec, (k + myW->getS() - 1) * myD, myD).vector;
    myW->AtWkV(&gv_sub, -k, R, &perm_col, myTempWkColRow);
    gsl_matrix_set_col(myDGammaTrMat, -k + myW->getS() - 1, &gv_sub);
  }

  int M = yr->size / myD;
  if (myD == 1) {
    dgamma_j_row = gsl_matrix_row(myDGammaTrMat, 0).vector;
    gsl_vector_add (myDGammaVec, &dgamma_j_row);
    tmv_prod_vector(myDGammaVec, myW->getS(), yr, M, res);  
  } else {
    gsl_vector_view res_stride = gsl_vector_subvector_with_stride(res, j, myD, M);       
    gsl_vector_view yr_stride = gsl_vector_subvector_with_stride(yr, j, myD, M);       
    gsl_matrix_view gamma_vec_mat = gsl_matrix_view_vector(myDGammaVec, 1, myDGammaVec->size);
    gsl_vector_memcpy(myTmpCol, &yr_stride.vector);
  
    gsl_vector_set_zero(res);
    tmv_prod_vector(myDGammaVec, myW->getS(), yr, M, &res_stride.vector);  
    tmv_prod_new(myDGammaTrMat, myW->getS(), myTmpCol, M, res, 1.0);  
  }
}

slraDGammaBBanded::slraDGammaBBanded( const slraSDependentStructure *s, int r ) :
   myD(s->getNplusD() - r), myW(s) {
}

slraDGammaBBanded::~slraDGammaBBanded(){
}

void slraDGammaBBanded::calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr ) {
                   
  for (size_t k = 0; k < myW->getM(); k++)  {
    
  }
}
