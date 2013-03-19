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

WLayeredHStructure::WLayeredHStructure( const double *m_l, size_t q, size_t n, 
    const double *w ) : myBase(m_l, q, n, NULL) {
  myInvWeights = gsl_vector_alloc(myBase.getNp());
  myInvSqrtWeights = gsl_vector_alloc(myBase.getNp());

  for (size_t l = 0; l < myInvWeights->size; l++) {
    if (!(w[l] > 0)) {
      throw new Exception("Value of weight not supported: %lf\n", w[l]);
    }
    gsl_vector_set(myInvWeights, l, (1 / w[l]));
    gsl_vector_set(myInvSqrtWeights, l, sqrt(gsl_vector_get(myInvWeights, l)));
  }
}

WLayeredHStructure::~WLayeredHStructure() {
  gsl_vector_free(myInvWeights);
  gsl_vector_free(myInvSqrtWeights);
}

void WLayeredHStructure::correctP( gsl_vector* p, gsl_matrix *R, 
                                   gsl_vector *yr, long wdeg ) {
  size_t l, k, sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix yr_matr = gsl_matrix_view_vector(yr, getN(), R->size2).matrix;
  gsl_vector yr_matr_row;
  gsl_vector_view res_sub, p_chunk_sub, inv_w_chunk;
  gsl_vector *res = gsl_vector_alloc(R->size1), 
             *weights = (wdeg == 2) ? myInvWeights : myInvSqrtWeights;

  for (k = 0; k < getN(); k++) {
    yr_matr_row = gsl_matrix_row(&yr_matr, k).vector; 
    gsl_blas_dgemv(CblasNoTrans, 1.0, R,  &yr_matr_row, 0.0, res); 
    
    for (l = 0, sum_np = 0, sum_nl = 0; l < getQ(); 
         sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
      res_sub = gsl_vector_subvector(res, sum_nl, getLayerLag(l));
      p_chunk_sub =  gsl_vector_subvector(p, k + sum_np, getLayerLag(l));
      if (wdeg != 0) {
        inv_w_chunk = gsl_vector_subvector(myInvWeights, k + sum_np, 
                                           getLayerLag(l));
        gsl_vector_mul(&res_sub.vector, &inv_w_chunk.vector);              
      }
      gsl_vector_sub(&p_chunk_sub.vector, &res_sub.vector); 
    }
  }

  gsl_vector_free(res);
}

void WLayeredHStructure::mulInvWij( gsl_matrix *matr, long i  ) const {
  size_t sum_np = i, sum_ml = 0, l, k;
  gsl_vector matr_row;

  for (l = 0; l < getQ(); sum_np += getLayerNp(l), sum_ml += getLayerLag(l), ++l) {
    for (k = 0; k < getLayerLag(l); ++k) {
      matr_row = gsl_matrix_row(matr, sum_ml + k).vector;
      gsl_vector_scale(&matr_row, getInvWeight(sum_np + k));
    }
  }
}

void WLayeredHStructure::WijB( gsl_matrix *res, long i, long j, 
         const gsl_matrix *B ) const {
  gsl_matrix_memcpy(res, B);
  if (i <= j) {
    mulInvWij(res, j);
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, 
        myBase.getWk(j-i), res);
  } else {
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1.0, 
        myBase.getWk(i-j), res);
    mulInvWij(res, i);
  }
}

void WLayeredHStructure::AtWijB( gsl_matrix *res, long i, long j, 
         const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *tmpWijB, 
         double beta ) const {
  gsl_matrix_scale(res, beta);
  size_t sum_np, ind_a, ind_b;
  size_t diff;
  
  diff = (j >= i ? j - i : i - j);
  ind_a = j - mymin(i, j);
  ind_b = i - mymin(i, j);

  for (size_t l = 0, sum_np = mymax(j, i); l < getQ(); 
       sum_np += getLayerNp(l), 
       ind_a += getLayerLag(l),
       ind_b += getLayerLag(l), ++l) {
    for (size_t k = 0; k + diff < getLayerLag(l); ++k) {
      const gsl_vector A_row = gsl_matrix_const_row(A, ind_a + k).vector;
      const gsl_vector B_row = gsl_matrix_const_row(B, ind_b + k).vector;
      gsl_blas_dger(getInvWeight(sum_np + k), &A_row, &B_row, res);
    }
  }
}   

void WLayeredHStructure::AtWijV( gsl_vector *res, long i, long j, 
         const gsl_matrix *A, const gsl_vector *V, 
                     gsl_vector *tmpWijV, double beta ) const {
  gsl_vector_scale(res, beta);
  size_t sum_np, ind_a, ind_v;
  size_t diff, l, k;

  diff = (j >= i ? j - i : i - j);
  ind_a = j - mymin(i, j);
  ind_v = i - mymin(i, j);
    
  for (l = 0, sum_np = mymax(j,i); l < getQ(); 
       sum_np += getLayerNp(l), 
       ind_a += getLayerLag(l),
       ind_v += getLayerLag(l), ++l) {
    for (k = 0; k + diff < getLayerLag(l); ++k) {
      const gsl_vector A_row = gsl_matrix_const_row(A, ind_a + k).vector;
      gsl_blas_daxpy(getInvWeight(sum_np + k) * 
                     gsl_vector_get(V, ind_v + k), &A_row, res);
    }
  }
}   

Cholesky *WLayeredHStructure::createCholesky( size_t D ) const {
  return new SDependentCholesky(this, D);
}

DGamma *WLayeredHStructure::createDGamma( size_t D ) const {
  return new SDependentDGamma(this, D);
}





