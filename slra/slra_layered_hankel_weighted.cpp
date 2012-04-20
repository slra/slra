#include <limits>
#include <memory.h>
#include <math.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

slraLayeredHankelWeightedStructure::
    slraLayeredHankelWeightedStructure( const double *oldNk, size_t q, int M, 
    const double *weights, bool block_weights ) : myBase(oldNk, q, M, NULL) {
  myInvSqrtWeights = gsl_vector_alloc(myBase.getNp());
  if (weights != NULL) {
    if (block_weights) {
      size_t sum_np = 0;
      for (size_t l = 0; l < myBase.getQ(); sum_np += myBase.getLayerNp(l), ++l) {
        if (!(weights[l] > 0)) {
          throw new slraException("This value of weight is not supported: %lf\n", weights[l]);
        }
        double inv_w = sqrt(1 / weights[l]);
        gsl_vector subw = 
            gsl_vector_subvector(myInvSqrtWeights, sum_np, myBase.getLayerNp(l)).vector;
        gsl_vector_set_all(&subw, inv_w);       
      }  
    } else {
      for (size_t l = 0; l < myInvSqrtWeights->size; l++) {
        if (!(weights[l] > 0)) {
          throw new slraException("This value of weight is not supported: %lf\n", weights[l]);
        }
        gsl_vector_set(myInvSqrtWeights, l, sqrt(1 / weights[l]));
      }
    }
  } else {
    gsl_vector_set_all(myInvSqrtWeights, 1);
  }
}

slraLayeredHankelWeightedStructure::~slraLayeredHankelWeightedStructure() {
  gsl_vector_free(myInvSqrtWeights);
}

void slraLayeredHankelWeightedStructure::
         correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  size_t l, k, sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix_view yr_matr = gsl_matrix_view_vector(yr, getM(), R->size2), b_xext;
  gsl_vector_view yr_matr_row, res_sub, p_chunk_sub, inv_w_chunk;
  gsl_vector *res = gsl_vector_alloc(R->size1);

  for (l = 0; l < getQ(); sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    b_xext = gsl_matrix_submatrix(R, sum_nl, 0, getLayerLag(l), R->size2); 
    res_sub = gsl_vector_subvector(res, 0, getLayerLag(l));
    for (k = 0; k < getM(); k++) {
      p_chunk_sub =  gsl_vector_subvector(p, k + sum_np, getLayerLag(l));
      inv_w_chunk = gsl_vector_subvector(myInvSqrtWeights, k + sum_np, getLayerLag(l));
      yr_matr_row = gsl_matrix_row(&yr_matr.matrix, k); 
      gsl_blas_dgemv(CblasNoTrans, 1.0, &b_xext.matrix, 
                      &yr_matr_row.vector, 0.0, &res_sub.vector); 
      gsl_vector_mul(&res_sub.vector, &inv_w_chunk.vector);              
      gsl_vector_sub(&p_chunk_sub.vector, &res_sub.vector); 
    }
  }
  gsl_vector_free(res);
}

void slraLayeredHankelWeightedStructure::mulInvWij( gsl_matrix *matr, int i  ) const {
  size_t sum_np = i, sum_nl = 0, l, k;
  double inv_w;
  gsl_vector matr_row;

  for (l = 0; l < getQ(); sum_np += getLayerNp(l), ++l) {
    for (k = 0; k < getLayerLag(l); ++k, ++sum_nl) {
      matr_row = gsl_matrix_row(matr, sum_nl).vector;
      inv_w = getInvSqrtWeights(sum_np + k);
      gsl_vector_scale(&matr_row, inv_w * inv_w);
    }
  }
}

void slraLayeredHankelWeightedStructure::WijB( gsl_matrix *res, int i, int j, 
         const gsl_matrix *B ) const {
  gsl_matrix_memcpy(res, B);
  if (j >= i) {
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, 
        myBase.getWk(j-i), res);
    mulInvWij(res, i);
  } else {
    mulInvWij(res, i);
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1.0, 
        myBase.getWk(i-j), res);
  }
}

void slraLayeredHankelWeightedStructure::AtWijB( gsl_matrix *res, int i, int j, 
         const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *tmpWijB, double beta ) const {
  WijB(tmpWijB, i, j, B);   
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, tmpWijB, beta, res);                   
}

void slraLayeredHankelWeightedStructure::AtWijV( gsl_vector *res, int i, int j,
                     const gsl_matrix *A, const gsl_vector *V, 
                     gsl_vector *tmpWijV, double beta ) const { /* Not optimal */
  gsl_matrix res_m = gsl_matrix_view_vector(res, res->size, 1).matrix;
  const gsl_matrix V_m = gsl_matrix_const_view_vector(V, V->size, 1).matrix;
  gsl_matrix tmpWijV_m = gsl_matrix_view_vector(tmpWijV, tmpWijV->size, 1).matrix;
  AtWijB(&res_m, i, j, A, &V_m, &tmpWijV_m, beta);
}

slraGammaCholesky *slraLayeredHankelWeightedStructure::
                       createGammaComputations( int r, double reg_gamma ) const {
  return new slraGammaCholeskyBBanded(this, r, reg_gamma);
}


slraDGamma *slraLayeredHankelWeightedStructure::
                createDerivativeComputations( int r ) const {
  return new slraDGammaBBanded(this, r);
}


