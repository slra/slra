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

HLayeredElWStructure::HLayeredElWStructure( const double *m_vec, size_t q, size_t n, 
    const double *w_vec ) : myBase(m_vec, q, n, NULL) {
  myInvWeights = gsl_vector_alloc(myBase.getNp());
  myInvSqrtWeights = gsl_vector_alloc(myBase.getNp());

  for (size_t l = 0; l < myInvWeights->size; l++) {
    if (!(w_vec[l] > 0)) {
      throw new Exception("Value of weight not supported: %lf\n", w_vec[l]);
    }
    gsl_vector_set(myInvWeights, l, (1 / w_vec[l]));
    gsl_vector_set(myInvSqrtWeights, l, sqrt(gsl_vector_get(myInvWeights, l)));
  }
}

HLayeredElWStructure::~HLayeredElWStructure() {
  gsl_vector_free(myInvWeights);
  gsl_vector_free(myInvSqrtWeights);
}

void HLayeredElWStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) {
  myBase.fillMatrixFromP(c, p);
}

void HLayeredElWStructure::multByGtUnweighted( gsl_vector* p, 
          const gsl_matrix *Rt, const gsl_vector *y, 
          double alpha, double beta, bool skipFixedBlocks ) {
  myBase.multByGtUnweighted(p, Rt, y, alpha, beta, skipFixedBlocks);        
}

void HLayeredElWStructure::multByWInv( gsl_vector* p, long deg ) const {
  if (deg == 0) {
    return;
  }
  if (deg == 1) {
    gsl_vector_mul(p, myInvSqrtWeights);              
  } else if (deg == 2) {
    gsl_vector_mul(p, myInvWeights);              
  }
}

void HLayeredElWStructure::mulInvWij( gsl_matrix *matr, long i_1  ) const {
  size_t sum_np = i_1, sum_ml = 0, l, k;
  gsl_vector matr_row;

  for (l = 0; l < getQ(); sum_np += getLayerNp(l), sum_ml += getLayerLag(l), ++l) {
    for (k = 0; k < getLayerLag(l); ++k) {
      matr_row = gsl_matrix_row(matr, sum_ml + k).vector;
      gsl_vector_scale(&matr_row, getInvWeight(sum_np + k));
    }
  }
}

void HLayeredElWStructure::VijB( gsl_matrix *X, long i_1, long j_1, 
         const gsl_matrix *B ) const {
  gsl_matrix_memcpy(X, B);
  if (i_1 <= j_1) {
    mulInvWij(X, j_1);
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, 
        myBase.getWk(j_1-i_1), X);
  } else {
    gsl_blas_dtrmm(CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1.0, 
        myBase.getWk(i_1-j_1), X);
    mulInvWij(X, i_1);
  }
}

void HLayeredElWStructure::AtVijB( gsl_matrix *X, long i_1, long j_1, 
         const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *tmpVijB, 
         double beta ) const {
  gsl_matrix_scale(X, beta);
  size_t sum_np, ind_a, ind_b;
  size_t diff;
  
  diff = (j_1 >= i_1 ? j_1 - i_1 : i_1 - j_1);
  ind_a = j_1 - mymin(i_1, j_1);
  ind_b = i_1 - mymin(i_1, j_1);

  for (size_t l = 0, sum_np = mymax(j_1, i_1); l < getQ(); 
       sum_np += getLayerNp(l), 
       ind_a += getLayerLag(l),
       ind_b += getLayerLag(l), ++l) {
    for (size_t k = 0; k + diff < getLayerLag(l); ++k) {
      const gsl_vector A_row = gsl_matrix_const_row(A, ind_a + k).vector;
      const gsl_vector B_row = gsl_matrix_const_row(B, ind_b + k).vector;
      gsl_blas_dger(getInvWeight(sum_np + k), &A_row, &B_row, X);
    }
  }
}   

void HLayeredElWStructure::AtVijV( gsl_vector *u, long i_1, long j_1, 
         const gsl_matrix *A, const gsl_vector *V, 
                     gsl_vector *tmpVijV, double beta ) const {
  gsl_vector_scale(u, beta);
  size_t sum_np, ind_a, ind_v;
  size_t diff, l, k;

  diff = (j_1 >= i_1 ? j_1 - i_1 : i_1 - j_1);
  ind_a = j_1 - mymin(i_1, j_1);
  ind_v = i_1 - mymin(i_1, j_1);
    
  for (l = 0, sum_np = mymax(j_1,i_1); l < getQ(); 
       sum_np += getLayerNp(l), 
       ind_a += getLayerLag(l),
       ind_v += getLayerLag(l), ++l) {
    for (k = 0; k + diff < getLayerLag(l); ++k) {
      const gsl_vector A_row = gsl_matrix_const_row(A, ind_a + k).vector;
      gsl_blas_daxpy(getInvWeight(sum_np + k) * 
                     gsl_vector_get(V, ind_v + k), &A_row, u);
    }
  }
}   

Cholesky *HLayeredElWStructure::createCholesky( size_t d ) const {
  return new MuDependentCholesky(this, d);
}

DGamma *HLayeredElWStructure::createDGamma( size_t d ) const {
  return new MuDependentDGamma(this, d);
}





