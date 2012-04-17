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

slraLayeredHankelWeightedStructure::slraLayeredHankelWeightedStructure( const double *oldNk, size_t q, int M, 
    const double *weights ) : slraLayeredHankelStructure(oldNk, q, M, NULL) {
  myInvSqrtWeights = gsl_vector_alloc(nvGetNp());
  if (weights != NULL) {
    for (size_t l = 0; l < myInvSqrtWeights->size; l++) {
      if (!(weights[l] > 0)) {
        throw new slraException("This value of weight is not supported: %lf\n", weights[l]);
      }
      gsl_vector_set(myInvSqrtWeights, l, sqrt(1 / weights[l]));
    }
  } else {
    gsl_vector_set_all(myInvSqrtWeights, 1);
  }
}

slraLayeredHankelWeightedStructure::~slraLayeredHankelWeightedStructure() {
  gsl_vector_free(myInvSqrtWeights);
}


void slraLayeredHankelWeightedStructure::correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  size_t l, k, sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix_view yr_matr = gsl_matrix_view_vector(yr, getM(), R->size2), b_xext;
  gsl_vector_view yr_matr_row, res_sub, p_chunk_sub, inv_w_chunk;
  gsl_vector *res = gsl_vector_alloc(R->size1);

  for (l = 0; l < getQ(); sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    b_xext = gsl_matrix_submatrix(R, sum_nl, 0, getLayerLag(l), R->size2); 
    res_sub = gsl_vector_subvector(res, 0, getLayerLag(l));
    if (!isLayerExact(l)) {    /* Subtract correction if needed */
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
  }
  
  gsl_vector_free(res);
}
