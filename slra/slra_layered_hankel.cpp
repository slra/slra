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

LayeredHStructure::LayeredHStructure( const double *oldNk, 
    size_t q, int M, const double *layer_w  ) :
                                      myQ(q), myM(M), mySA(NULL)  {
  mySA = new Layer[myQ];
  for (size_t l = 0; l < myQ; l++) {
    mySA[l].blocks_in_row = oldNk[l];
    if (layer_w != NULL && !(layer_w[l] > 0)) {
      throw new Exception("This value of weight is not supported: %lf\n", layer_w[l]);
    }
    mySA[l].inv_w = (layer_w != NULL) ? (1 / layer_w[l]) : 1;
  }    
   
  computeStats(); 
  computeWkParams();
}

LayeredHStructure::~LayeredHStructure() {
  if (mySA != NULL) {
    delete[] mySA;
  }
  if (myA != NULL) {
    for (int k = 0; k < myMaxLag; k++) {
      gsl_matrix_free(myA[k]);
    }
    free(myA);
  }
}

void LayeredHStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  {
  size_t sum_np = 0, sum_nl = 0, l, j;
  gsl_matrix_view c_chunk, c_chunk_sub;
 
  for (l = 0; l < getQ(); sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    c_chunk = gsl_matrix_submatrix(c, 0, sum_nl, getM(), getLayerLag(l));
    for (j = 0; j < getLayerLag(l); j++) {
      gsl_vector_const_view psub = gsl_vector_const_subvector(p, sum_np + j, getM());
      gsl_matrix_set_col(&c_chunk.matrix, j, &psub.vector);
    }  
  }
}

void LayeredHStructure::computeWkParams() {
  int k, l, i, imax, sum_nl, rep;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;

  myA = (gsl_matrix**) malloc(getMaxLag() * sizeof(gsl_matrix *));
  /* construct w */
  for (k = 0; k < getMaxLag(); k++) { 
    zk   = gsl_matrix_alloc(getNplusD(), getNplusD());
    gsl_matrix_set_zero(zk);

    for (sum_nl = 0, l = 0; l < getQ(); sum_nl += getLayerLag(l), ++l) { 
      zkl = gsl_matrix_submatrix(zk, sum_nl, sum_nl, getLayerLag(l), getLayerLag(l)); 
      if (k < getLayerLag(l)) {
        gsl_vector_view diag = gsl_matrix_subdiagonal(&zkl.matrix, k);
        gsl_vector_set_all(&diag.vector, getLayerInvWeight(l)); 
      }
    }
    myA[k] = zk;
  }
}

void LayeredHStructure::computeStats() {
  int l;
  for (l = 0, myNplusD = 0, myMaxLag = 1; l < myQ; myNplusD += getLayerLag(l), ++l) {
    if ((!isLayerExact(l)) && getLayerLag(l) > myMaxLag) {
      myMaxLag = mySA[l].blocks_in_row;
    }
  }
}

void LayeredHStructure::correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  size_t l, k, sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix_view yr_matr = gsl_matrix_view_vector(yr, getM(), R->size2), b_xext;
  gsl_vector_view yr_matr_row, res_sub, p_chunk_sub;
  gsl_vector *res = gsl_vector_alloc(R->size1);

  for (l = 0; l < getQ(); sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    b_xext = gsl_matrix_submatrix(R, sum_nl, 0, getLayerLag(l), R->size2); 
    res_sub = gsl_vector_subvector(res, 0, getLayerLag(l));
    if (!isLayerExact(l)) {    /* Subtract correction if needed */
      for (k = 0; k < getM(); k++) {
        p_chunk_sub =  gsl_vector_subvector(p, k + sum_np, getLayerLag(l));
        yr_matr_row = gsl_matrix_row(&yr_matr.matrix, k); 
        gsl_blas_dgemv(CblasNoTrans, sqrt(getLayerInvWeight(l)), &b_xext.matrix, 
                        &yr_matr_row.vector, 0.0, &res_sub.vector); 
        gsl_vector_sub(&p_chunk_sub.vector, &res_sub.vector); 
      }
    }
  }
  
  gsl_vector_free(res);
}

Cholesky *LayeredHStructure::createCholesky( int D, double reg_gamma ) const {
#ifdef USE_SLICOT 
  return new StationaryCholeskySlicot(this, D, reg_gamma);
#else  /* USE_SLICOT */
  return new StationaryCholesky(this, D, reg_gamma);
#endif /* USE_SLICOT */
}

DGamma *LayeredHStructure::createDGamma( int D ) const {
  return new StationaryDGamma(this, D);
}


void LayeredHStructure::setM( int m ) {
  myM = m <= 0 ? 1 : m;
}

void LayeredHStructure::WkB( gsl_matrix *res, int k, const gsl_matrix *B ) const {

  gsl_matrix_memcpy(res, B);
  gsl_blas_dtrmm(CblasLeft, CblasLower, (k > 0 ? CblasNoTrans : CblasTrans), CblasNonUnit, 1.0, getWk(abs(k)), res);
}

void LayeredHStructure::AtWkB( gsl_matrix *res, int k,
         const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *tmpWkB, double beta ) const {
  WkB(tmpWkB, k, B);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, tmpWkB, beta, res);       
}

void LayeredHStructure::AtWkV( gsl_vector *res, int k, 
         const gsl_matrix *A, const gsl_vector *V, gsl_vector *tmpWkV, double beta ) const {
  gsl_blas_dcopy(V, tmpWkV);
  gsl_blas_dtrmv(CblasLower, (k > 0 ? CblasNoTrans : CblasTrans), CblasNonUnit, getWk(abs(k)), tmpWkV);
  gsl_blas_dgemv(CblasTrans, 1.0, A, tmpWkV, beta, res);       
}

Cholesky *MosaicHStructure::createCholesky( int D, double reg_gamma ) const {
  if (myWkIsCol) {
    return new SameStripedStationaryCholesky(this, D, 1, reg_gamma);
  } else {
    return new StripedCholesky(this, D, reg_gamma);
  }
}

typedef Structure* pStructure;

pStructure * MosaicHStructure::allocStripe(gsl_vector *oldNk, gsl_vector *oldMl,  
                gsl_vector *Wk,  bool wkIsCol )  {
  pStructure *res = new pStructure[oldMl->size];

  for (size_t k = 0; k < oldMl->size; k++) {
    res[k] = new LayeredHStructure(oldNk->data, oldNk->size, oldMl->data[k], 
      (Wk != NULL ? Wk->data : NULL) );
    if (Wk != NULL && wkIsCol) {
      Wk->data += oldNk->size;
    }
  }

  return res;
}

MosaicHStructure::MosaicHStructure( gsl_vector *oldNk, gsl_vector *oldMl,  
          gsl_vector *Wk, bool wkIsCol ) :
        StripedStructure(oldMl->size, allocStripe(oldNk, oldMl, Wk)),
        myWkIsCol(wkIsCol)  {
}









