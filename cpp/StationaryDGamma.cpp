#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

StationaryDGamma::StationaryDGamma( const StationaryStructure *s, size_t D ) :
    myD(D), myW(s) {
  myTempVkColRow = gsl_vector_alloc(myW->getM());
  myDGammaVec = gsl_vector_alloc(myD * (2 * myW->getMu() - 1));
  myDGammaTrMat = gsl_matrix_alloc(myD, 2 * myW->getMu() - 1);
  myDGamma = gsl_matrix_alloc(myD, myD * (2 * myW->getMu() - 1));
  myTmpCol = gsl_vector_alloc(myW->getN());
  myVk_R =  gsl_matrix_alloc(myW->getM(), myD);
  myN_k = gsl_matrix_alloc(myD, myD);
  myEye = gsl_matrix_alloc(myW->getM(), myW->getM());
  gsl_matrix_set_identity(myEye);
}

StationaryDGamma::~StationaryDGamma() {
  gsl_vector_free(myTempVkColRow);
  gsl_vector_free(myTmpCol);
  gsl_vector_free(myDGammaVec);
  gsl_matrix_free(myDGammaTrMat);
  gsl_matrix_free(myDGamma);
  gsl_matrix_free(myVk_R);
  gsl_matrix_free(myN_k);
  gsl_matrix_free(myEye);
}

void StationaryDGamma::calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt, 
                                      const gsl_matrix *Yt ) {
  size_t n = Yt->size1;
  gsl_matrix YrB, YrT;

  gsl_matrix_set_zero(At);
  
  for (size_t k = 0; k < myW->getMu(); k++) {
    YrT = gsl_matrix_const_submatrix(Yt, 0, 0, n - k, myD).matrix;
    YrB = gsl_matrix_const_submatrix(Yt, k, 0, n - k, myD).matrix;
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &YrB, &YrT, 0.0, myN_k);

    myW->VkB(myVk_R, k, Rt);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 2.0, myVk_R, myN_k, 1.0, At);

    if (k > 0) {
      myW->VkB(myVk_R, -k, Rt);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 2.0, myVk_R, myN_k, 1.0, At);
    }
  }    
}

void StationaryDGamma::calcDijGammaYr( gsl_vector *z, const gsl_matrix *R, 
         size_t j_1, size_t i_1, const gsl_vector *y ) {
  gsl_vector gv_sub, e_j = gsl_matrix_column(myEye, j_1).vector, dgammajrow,
             res_stride, y_stride;
  long S = myW->getMu();           

  for (long k = 1 - S; k < S; k++) {
    gv_sub = gsl_vector_subvector(myDGammaVec, (k + S - 1) * myD, 
              myD).vector;
    myW->AtVkV(&gv_sub, -k, R, &e_j, myTempVkColRow);
    gsl_matrix_set_col(myDGammaTrMat, -k + myW->getMu() - 1, &gv_sub);
  }

  size_t n = y->size / myD;
  if (myD == 1) {
    dgammajrow = gsl_matrix_row(myDGammaTrMat, 0).vector;
    gsl_vector_add (myDGammaVec, &dgammajrow);
    tmv_prod_vector(myDGammaVec, myW->getMu(), y, n, z);  
  } else {
    res_stride = gsl_vector_subvector_with_stride(z, i_1, myD, n).vector;       
    y_stride = gsl_vector_const_subvector_with_stride(y, i_1, myD, n).vector;       
    gsl_matrix_view gamma_vec_mat = gsl_matrix_view_vector(myDGammaVec, 1, 
                                                           myDGammaVec->size);
    gsl_vector_memcpy(myTmpCol, &y_stride);
  
    gsl_vector_set_zero(z);
    tmv_prod_vector(myDGammaVec, myW->getMu(), y, n, &res_stride);  
    tmv_prod_new(myDGammaTrMat, myW->getMu(), myTmpCol, n, z, 1.0);  
  }
}


