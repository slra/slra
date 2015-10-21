#include <memory.h>
#include "slra.h"

void tmv_prod_vector( gsl_vector *T, size_t s, const gsl_vector* v, size_t m,
                     gsl_vector* p ) {
  double res;
  size_t i, s_1 = s - 1;
  size_t bCols = T->size / (2 * s - 1);
  size_t j1 = 0, t1 = s_1 * bCols, len = GSL_MIN(s, m) * bCols;
  gsl_vector subT, subv; /* subvectors of v and p */
  
  for (i = 0; i < m; ++i) {
    subv = gsl_vector_const_subvector(v, j1, len).vector;
    subT = gsl_vector_subvector(T, t1, len).vector;
    
    gsl_blas_ddot(&subT, &subv, &res);
    gsl_vector_set(p, i, res);
    
    if (i < s_1) {
      t1 -= bCols;
    } else {
      j1 += bCols;
      len -= bCols;
    }
    if (i + s < m) {
      len += bCols;
    }
  }
}

void tmv_prod_new( gsl_matrix *T, size_t s, const gsl_vector* v,
                   size_t m, gsl_vector* p, double beta ) {
  size_t i, s_1 = s - 1;
  size_t D = T->size1;
  size_t bCols = T->size2 / (2 * s - 1);
  size_t j1 = 0, t1 = s_1 * bCols, len = GSL_MIN(s, m) * bCols;
  gsl_matrix subT;
  gsl_vector subv, subp; 	/* subvectors of v and p */

  for (i = 0; i < m; ++i) {
    
    subp = gsl_vector_subvector(p, i * D, D).vector;
    subv = gsl_vector_const_subvector(v, j1, len).vector;
    subT = gsl_matrix_submatrix(T, 0, t1, D, len).matrix;
    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &subT, &subv, beta, &subp);
    
    if (i < s_1) {
      t1 -= bCols;
    } else {
      j1 += bCols;
      len -= bCols;
    }
    if (i + s < m) {
      len += bCols;
    }
  }
}

void copyLowerTrg( gsl_matrix * dest, const gsl_matrix *src  ) {
  for (size_t i = 0; i < dest->size1; i++) {
    for (size_t j = 0; (j < dest->size2) && (j <= i); j++) {
      gsl_matrix_set(dest, i, j, gsl_matrix_get(src, i, j));
    }
  }
}


void shiftLowerTrg( gsl_matrix * dest, const gsl_matrix *src  ) {
  gsl_matrix_set_zero(dest);
  for (size_t i = 0; i < dest->size1; i++) {
    for (size_t j = 0; (j < dest->size2) && (j <= i); j++) {
      gsl_matrix_set(dest, i, (dest->size2 - 1 - i) + j, 
          gsl_matrix_get(src, i, j));
    }
  }
}

void id_kron_a( const gsl_matrix *A, size_t d,  gsl_matrix *IkronA ) {
  gsl_matrix_set_zero(IkronA);
  for (size_t  k = 0; k < d; k++) {
    gsl_matrix subA = gsl_matrix_submatrix(IkronA, k * A->size1, k * A->size2, 
                           A->size1, A->size2).matrix;
    gsl_matrix_memcpy(&subA, A);
  }
}

void ls_solve( const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *X ) {
  size_t d = B->size2, i, j, one = 1, lwork = -1, info;
  gsl_matrix *IkronA = gsl_matrix_alloc(A->size1 * d, A->size2 * d);
  size_t vecBsize = B->size2 * B->size1; 
  double tmp;
  
  double *vecIkronA = new double[IkronA->size1 * IkronA->size2],
         *vecB = new double[B->size1 * B->size2]; 
                          
  id_kron_a(A, d, IkronA);                                      
  gsl_matrix_vectorize(vecIkronA, IkronA);
  gsl_matrix_vectorize(vecB, B);
  
  dgels_("N", &IkronA->size1, &IkronA->size2, &one, vecIkronA,
         &IkronA->size1, vecB, &vecBsize, &tmp, &lwork, &info);
  double *work = new double[lwork = tmp];
  dgels_("N", &IkronA->size1, &IkronA->size2, &one, vecIkronA,
         &IkronA->size1, vecB, &vecBsize, work, &lwork, &info);

  delete [] work;
  
  gsl_matrix_free(IkronA);
  gsl_matrix_vec_inv(X, vecB);
  
  delete [] vecB;
  delete [] vecIkronA;
}
