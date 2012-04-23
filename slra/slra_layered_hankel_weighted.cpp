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

WLayeredHStructure::WLayeredHStructure( const double *oldNk, size_t q, int M, 
    const double *weights ) : myBase(oldNk, q, M, NULL) {
    
  myInvSqrtWeights = gsl_vector_alloc(myBase.getNp());
  if (weights != NULL) {
    for (size_t l = 0; l < myInvSqrtWeights->size; l++) {
      if (!(weights[l] > 0)) {
        throw new Exception("Value of weight not supported: %lf\n", weights[l]);
      }
      gsl_vector_set(myInvSqrtWeights, l, sqrt(1 / weights[l]));
    }
  } else {
    gsl_vector_set_all(myInvSqrtWeights, 1);
  }
}

WLayeredHStructure::~WLayeredHStructure() {
  gsl_vector_free(myInvSqrtWeights);
}

void WLayeredHStructure::correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  size_t l, k, sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix_view yr_matr = gsl_matrix_view_vector(yr, getM(), R->size2), b_xext;
  gsl_vector_view yr_matr_row, res_sub, p_chunk_sub, inv_w_chunk;
  gsl_vector *res = gsl_vector_alloc(R->size1);

  for (l = 0; l < getQ(); sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    b_xext = gsl_matrix_submatrix(R, sum_nl, 0, getLayerLag(l), R->size2); 
    res_sub = gsl_vector_subvector(res, 0, getLayerLag(l));
    for (k = 0; k < getM(); k++) {
      p_chunk_sub =  gsl_vector_subvector(p, k + sum_np, getLayerLag(l));
      inv_w_chunk = gsl_vector_subvector(myInvSqrtWeights, k + sum_np, 
                                         getLayerLag(l));

//      gsl_vector_memcpy(&p_chunk_sub.vector, &inv_w_chunk.vector);
//      print_arr(inv_w_chunk.vector.data, inv_w_chunk.vector.size);
     
      yr_matr_row = gsl_matrix_row(&yr_matr.matrix, k); 
      gsl_blas_dgemv(CblasNoTrans, 1.0, &b_xext.matrix, 
                      &yr_matr_row.vector, 0.0, &res_sub.vector); 
      gsl_vector_mul(&res_sub.vector, &inv_w_chunk.vector);              
      gsl_vector_sub(&p_chunk_sub.vector, &res_sub.vector); 
    }
  }
  
//  gsl_vector_memcpy(p, myInvSqrtWeights);
//  PRINTF("InvWeights\n");
//  print_vec(myInvSqrtWeights);

  gsl_vector_free(res);
}

void WLayeredHStructure::mulInvWij( gsl_matrix *matr, int i  ) const {
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

void WLayeredHStructure::WijB( gsl_matrix *res, int i, int j, 
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

void WLayeredHStructure::AtWijB( gsl_matrix *res, int i, int j, 
         const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *tmpWijB, double beta ) const {
  gsl_matrix_scale(res, beta);
  size_t sum_np, sum_nl = 0;
  int l, k;
  double inv_w;

  if (i > j) {
    int tmp;
    const gsl_matrix* C;
    tmp = i; i = j; j = tmp;

    C = A; A = B; B = C;  
  }

  int diff = j - i;
  sum_nl = 0;
  for (l = 0, sum_np = j; l < getQ(); 
       sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    for (k = 0; k < ((int)getLayerLag(l)) - diff; ++k) {
      const gsl_vector A_row = gsl_matrix_const_row(A, sum_nl + k + diff).vector;
      const gsl_vector B_row = gsl_matrix_const_row(B, sum_nl + k).vector;
      inv_w = getInvSqrtWeights(sum_np + k);
      gsl_blas_dger(inv_w * inv_w, &A_row, &B_row, res);
    }
  }
}   

void WLayeredHStructure::AtWijV( gsl_vector *res, int i, int j, 
         const gsl_matrix *A, const gsl_vector *V, 
                     gsl_vector *tmpWijV, double beta ) const {
  gsl_vector_scale(res, beta);
  
  size_t sum_np, sum_nl = 0;
  int l, k;
  double inv_w;

  if (j >= i) {
    int diff = j - i;
    sum_nl = 0;
    for (l = 0, sum_np = j; l < getQ(); 
         sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
      for (k = 0; k < ((int)getLayerLag(l)) - diff; ++k) {
        const gsl_vector A_row = gsl_matrix_const_row(A, sum_nl + k + diff).vector;
        inv_w = getInvSqrtWeights(sum_np + k);
        gsl_blas_daxpy(inv_w * inv_w * gsl_vector_get(V, sum_nl + k), &A_row, res);
      }
    }
  } else {
    sum_nl = 0;
    int diff = i - j; 
    for (l = 0, sum_np = i; l < getQ(); 
         sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
      for (k = 0; k < ((int)getLayerLag(l)) - diff; ++k) {
        const gsl_vector A_row = gsl_matrix_const_row(A, sum_nl + k).vector;
        inv_w = getInvSqrtWeights(sum_np + k);
        gsl_blas_daxpy(inv_w * inv_w * gsl_vector_get(V, sum_nl + k + diff), &A_row, res);
      }
    }
  }
}   

Cholesky *WLayeredHStructure::createCholesky( int D, double reg_gamma ) const {
  return new SDependentCholesky(this, D, reg_gamma);
}


DGamma *WLayeredHStructure::createDGamma( int D ) const {
  return new SDependentDGamma(this, D);
}


typedef Structure* pStructure;

pStructure *WMosaicHStructure::allocStripe( size_t q, size_t N, 
     double *oldNk,  double *oldMl, double *Wk )  {
  pStructure *res = new pStructure[N];

  for (size_t k = 0; k < N; k++) {
    res[k] = new WLayeredHStructure(oldNk, q, oldMl[k], Wk);
    if (Wk != NULL) {
//      print_arr(Wk, res[k]->getNp());
      Wk += res[k]->getNp();
    }
  }

  return res;
}

WMosaicHStructure::WMosaicHStructure( size_t q, size_t N,
    double *oldNk, double *oldMl, double *Wk ) :
        StripedStructure(N, allocStripe(q, N, oldNk, oldMl, Wk)) {
}


