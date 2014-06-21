#include <limits>
#include <memory.h>
#include <cstdarg>
#include "slra.h"

HLayeredBlWStructure::HLayeredBlWStructure( const double *m_vec, 
    size_t q, size_t n, const double *w_vec  ) : myQ(q), myN(n), mySA(NULL)  {
  mySA = new Layer[myQ];
 
  for (size_t l_1 = 0; l_1 < myQ; ++l_1) {
    mySA[l_1].blocks_in_row = m_vec[l_1];
    if (w_vec != NULL && !(w_vec[l_1] > 0)) {
      throw new Exception("This value of weight is not supported: %lf\n", w_vec[l_1]);
    }
    mySA[l_1].inv_w = (w_vec != NULL) ? (1 / w_vec[l_1]) : 1.0;
  }    
   
  computeStats(); 
  computeVkParams();
}

HLayeredBlWStructure::~HLayeredBlWStructure() {
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

void HLayeredBlWStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) {
  size_t sum_np = 0, sum_nl = 0, l_1, j;
  gsl_vector psub;
 
  for (l_1 = 0; l_1 < getQ(); sum_np += getLayerNp(l_1), 
                              sum_nl += getLayerLag(l_1), ++l_1) {
    for (j = 0; j < getLayerLag(l_1); ++j) {
      psub = gsl_vector_const_subvector(p, sum_np + j, getN()).vector;
      gsl_matrix_set_col(c, j + sum_nl, &psub);
    }  
  }
}

void HLayeredBlWStructure::computeVkParams() {
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

void HLayeredBlWStructure::computeStats() {
  size_t l;
  for (l = 0, myM = 0, myMaxLag = 1; l < myQ; 
       myM += getLayerLag(l), ++l) {
    if ((!isLayerExact(l)) && getLayerLag(l) > myMaxLag) {
      myMaxLag = getLayerLag(l);
    }
  }
}

void HLayeredBlWStructure::multByGtUnweighted( gsl_vector* p, 
          const gsl_matrix *Rt, const gsl_vector *y, 
          double alpha, double beta, bool skipFixedBlocks ) {
  size_t l, k, sum_np = 0, sum_nl = 0, D = Rt->size2;
  gsl_matrix Y = gsl_matrix_const_view_vector(y,getN(), D).matrix, RtSub;
  gsl_vector Y_row, psub;

  for (l = 0; l < getQ(); 
       sum_np += getLayerNp(l), sum_nl += getLayerLag(l), ++l) {
    RtSub = gsl_matrix_const_submatrix(Rt, sum_nl, 0, getLayerLag(l), D).matrix; 
              
    if (!(skipFixedBlocks && isLayerExact(l))) {  
      for (k = 0; k < getN(); k++) {
        psub = gsl_vector_subvector(p, k + sum_np, getLayerLag(l)).vector;
        Y_row = gsl_matrix_row(&Y, k).vector; 
        gsl_blas_dgemv(CblasNoTrans, alpha, &RtSub, &Y_row, beta, &psub); 
      }
    } 
  }
} 
void HLayeredBlWStructure::multByWInv( gsl_vector* p, long deg ) {
  size_t l_1, k, sum_np = 0;
  gsl_vector psub;
  
  if (deg == 0) {
    return;
  }
  for (l_1 = 0; l_1 < getQ(); sum_np += getLayerNp(l_1), ++l_1) {
    psub = gsl_vector_subvector(p, sum_np, getLayerNp(l_1)).vector;
    gsl_vector_scale(&psub, (deg == 2) ? getLayerInvWeight(l_1):
                                         sqrt(getLayerInvWeight(l_1)));
  }
}

Cholesky *HLayeredBlWStructure::createCholesky( size_t d ) const {
#ifdef USE_SLICOT 
  return new StationaryCholeskySlicot(this, d);
#else  /* USE_SLICOT */
  return new StationaryCholesky(this, d);
#endif /* USE_SLICOT */
}

DGamma *HLayeredBlWStructure::createDGamma( size_t d ) const {
  return new StationaryDGamma(this, d);
}

void HLayeredBlWStructure::VkB( gsl_matrix *X, long k, 
                                const gsl_matrix *B ) const {
  gsl_matrix_memcpy(X, B);
  gsl_blas_dtrmm(CblasLeft, CblasLower, (k > 0 ? CblasNoTrans : CblasTrans), 
                 CblasNonUnit, 1.0, getWk(abs(k)), X);
}

void HLayeredBlWStructure::AtVkB( gsl_matrix *X, long k, const gsl_matrix *A,
         const gsl_matrix *B, gsl_matrix *tmpVkB, double beta ) const {
  VkB(tmpVkB, k, B);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, tmpVkB, beta, X);       
}

void HLayeredBlWStructure::AtVkV( gsl_vector *u, long k, const gsl_matrix *A,
         const gsl_vector *v, gsl_vector *tmpVkV, double beta ) const {
  gsl_blas_dcopy(v, tmpVkV);
  gsl_blas_dtrmv(CblasLower, (k > 0 ? CblasNoTrans : CblasTrans),
                 CblasNonUnit, getWk(abs(k)), tmpVkV);
  gsl_blas_dgemv(CblasTrans, 1.0, A, tmpVkV, beta, u);       
}






