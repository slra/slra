#include <limits>
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

LayeredHStructure::LayeredHStructure( const double *m_l, 
    size_t q, size_t n, const double *w_l  ) : myQ(q), myN(n), mySA(NULL)  {
  mySA = new Layer[myQ];
 
  for (size_t l = 0; l < myQ; l++) {
    mySA[l].blocks_in_row = m_l[l];
    if (w_l != NULL && !(w_l[l] > 0)) {
      throw new Exception("This value of weight is not supported: %lf\n", 
                          w_l[l]);
    }
    mySA[l].inv_w = (w_l != NULL) ? (1 / w_l[l]) : 1.0;
  }    
   
  computeStats(); 
  computeWkParams();
}

LayeredHStructure::~LayeredHStructure() {
  if (mySA != NULL) {
    delete[] mySA;
  }
  if (myA != NULL) {
    for (size_t k = 0; k < myMaxLag; k++) {
      gsl_matrix_free(myA[k]);
    }
    free(myA);
  }
}

void LayeredHStructure::fillMatrixFromP( gsl_matrix* c, 
                                         const gsl_vector* p ) {
  size_t sum_np = 0, sum_nl = 0, l, j;
  gsl_matrix_view c_chunk, c_chunk_sub;
 
  for (l = 0; l < getQ(); 
       sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    c_chunk = gsl_matrix_submatrix(c, 0, sum_nl, getN(), getLayerLag(l));
    for (j = 0; j < getLayerLag(l); j++) {
      gsl_vector_const_view psub = gsl_vector_const_subvector(p, sum_np + j, 
                                                              getN());
      gsl_matrix_set_col(&c_chunk.matrix, j, &psub.vector);
    }  
  }
}

void LayeredHStructure::computeWkParams() {
  size_t k, l, i, imax, sum_nl, rep;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;

  myA = (gsl_matrix**) malloc(getMaxLag() * sizeof(gsl_matrix *));
  
  /* construct w */
  for (k = 0; k < getMaxLag(); k++) { 
    zk   = gsl_matrix_alloc(getM(), getM());
    gsl_matrix_set_zero(zk);

    for (sum_nl = 0, l = 0; l < getQ(); sum_nl += getLayerLag(l), ++l) { 
      zkl = gsl_matrix_submatrix(zk, sum_nl, sum_nl, getLayerLag(l), 
                getLayerLag(l)); 
      if (k < getLayerLag(l)) {
        gsl_vector_view diag = gsl_matrix_subdiagonal(&zkl.matrix, k);
        gsl_vector_set_all(&diag.vector, getLayerInvWeight(l)); 
      }
    }
    myA[k] = zk;
  }
}

void LayeredHStructure::computeStats() {
  size_t l;
  for (l = 0, myM = 0, myMaxLag = 1; l < myQ; 
       myM += getLayerLag(l), ++l) {
    if ((!isLayerExact(l)) && getLayerLag(l) > myMaxLag) {
      myMaxLag = getLayerLag(l);
    }
  }
}

void LayeredHStructure::correctP( gsl_vector* p, gsl_matrix *R, 
                                  gsl_vector *yr, int wdeg ) {
  size_t l, k, sum_np = 0, sum_nl = 0, p_len, D = R->size2;
  gsl_matrix yr_matr = gsl_matrix_view_vector(yr, getN(), D).matrix, b_xext;
  gsl_vector yr_matr_row, res_sub, p_chunk_sub;
  gsl_vector *res = gsl_vector_alloc(R->size1);
  double w_scale;

  for (l = 0; l < getQ(); 
       sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    b_xext = gsl_matrix_submatrix(R, sum_nl, 0, getLayerLag(l), D).matrix; 
    res_sub = gsl_vector_subvector(res, 0, getLayerLag(l)).vector;
    w_scale = (wdeg == 2) ? getLayerInvWeight(l) :
                (wdeg == 1 ? sqrt(getLayerInvWeight(l)) : 1.0);
    if (!isLayerExact(l)) {    /* Subtract correction if needed */
      for (k = 0; k < getN(); k++) {
        p_chunk_sub = gsl_vector_subvector(p, k + sum_np, 
                                           getLayerLag(l)).vector;
        yr_matr_row = gsl_matrix_row(&yr_matr, k).vector; 
        gsl_blas_dgemv(CblasNoTrans, w_scale, &b_xext, &yr_matr_row, 0.0,
            &res_sub); 
        gsl_vector_sub(&p_chunk_sub, &res_sub); 
      }
    }
  }
  
  gsl_vector_free(res);
}

Cholesky *LayeredHStructure::createCholesky( size_t D ) const {
#ifdef USE_SLICOT 
  return new StationaryCholeskySlicot(this, D);
#else  /* USE_SLICOT */
  return new StationaryCholesky(this, D);
#endif /* USE_SLICOT */
}

DGamma *LayeredHStructure::createDGamma( size_t D ) const {
  return new StationaryDGamma(this, D);
}


void LayeredHStructure::WkB( gsl_matrix *res, int k, 
                             const gsl_matrix *B ) const {
  gsl_matrix_memcpy(res, B);
  gsl_blas_dtrmm(CblasLeft, CblasLower, (k > 0 ? CblasNoTrans : CblasTrans), 
                 CblasNonUnit, 1.0, getWk(abs(k)), res);
}

void LayeredHStructure::AtWkB( gsl_matrix *res, int k, const gsl_matrix *A,
         const gsl_matrix *B, gsl_matrix *tmpWkB, double beta ) const {
  WkB(tmpWkB, k, B);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, tmpWkB, beta, res);       
}

void LayeredHStructure::AtWkV( gsl_vector *res, int k, const gsl_matrix *A,
         const gsl_vector *V, gsl_vector *tmpWkV, double beta ) const {
  gsl_blas_dcopy(V, tmpWkV);
  gsl_blas_dtrmv(CblasLower, (k > 0 ? CblasNoTrans : CblasTrans),
                 CblasNonUnit, getWk(abs(k)), tmpWkV);
  gsl_blas_dgemv(CblasTrans, 1.0, A, tmpWkV, beta, res);       
}

Cholesky *MosaicHStructure::createCholesky( size_t D ) const {
  if (myWkIsCol) {
    return new SameStripedStationaryCholesky(this, D, 1);
  } else {
    return new StripedCholesky(this, D);
  }
}

typedef Structure* pStructure;

pStructure * MosaicHStructure::allocStripe( gsl_vector *m_l,
                 gsl_vector *n_k, gsl_vector *w )  {
  pStructure *res = new pStructure[n_k->size];

  for (size_t k = 0; k < n_k->size; k++) {
    res[k] = new LayeredHStructure(m_l->data, m_l->size, n_k->data[k], 
                                   (w != NULL ? w->data : NULL));
    if (w != NULL && (w->size != m_l->size)) {
      w->data += m_l->size;
    }
  }
  return res;
}

MosaicHStructure::MosaicHStructure( gsl_vector *m_l, gsl_vector *n_k,  
     gsl_vector *w ) : StripedStructure(n_k->size, allocStripe(m_l, n_k, w)),
        myWkIsCol(w == NULL || w->size == m_l->size)  {
}









