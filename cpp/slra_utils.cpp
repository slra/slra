#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

void tmv_prod_vector( gsl_vector *T, size_t s, const gsl_vector* v, size_t m, 
                      gsl_vector* p ) {
  double res;
  size_t i, temp, s_1 = s - 1;
  size_t bCols = T->size / (2 * s - 1), corner_row_lim = GSL_MIN(s_1, m/2); 
  gsl_vector subT, subv; /* subvectors of v and p */

  for (i = 0; i < corner_row_lim; i++) {  /* beginning and end parts */
    temp = GSL_MIN(s + i, m) * bCols;
    /* beginning part of the product p */
    subv = gsl_vector_const_subvector(v, 0, temp).vector;
    subT = gsl_vector_subvector(T, (s_1 - i) * bCols, temp).vector;
    gsl_blas_ddot(&subT, &subv, &res);
    gsl_vector_set(p, i, res);
    
    /* last part of the product p */
    subv = gsl_vector_const_subvector(v, v->size - temp, temp).vector;
    subT = gsl_vector_subvector(T, (s + i) * bCols - temp, temp).vector;    
    gsl_blas_ddot(&subT, &subv, &res);
    gsl_vector_set(p, p->size - (i + 1), res);
  }
  /* middle part */
  for (i = s_1; i < m - s_1; i++) {
    subv = gsl_vector_const_subvector(v, (i - s_1) * bCols, T->size).vector;
    gsl_blas_ddot(T, &subv, &res);
    gsl_vector_set(p, i,  res);
  }
}

void tmv_prod_new( gsl_matrix *T, size_t s, const gsl_vector* v, 
                   size_t m, gsl_vector* p, double beta ) {
  size_t i, temp, s_1 = s - 1;
  size_t D = T->size1;
  size_t bCols = T->size2 / (2 * s - 1), corner_row_lim = GSL_MIN(s_1, m/2); 
  gsl_matrix subT;
  gsl_vector subv, subp; 	/* subvectors of v and p */

  for (i = 0; i < corner_row_lim; i++) {  /* beginning and end parts */
    temp = GSL_MIN(s + i, m) * bCols;
    /* beginning part of the product p */
    subp = gsl_vector_subvector(p, i * D, D).vector;
    subv = gsl_vector_const_subvector(v, 0, temp).vector;
    subT = gsl_matrix_submatrix(T, 0, (s_1 - i) * bCols, D, temp).matrix;
    gsl_blas_dgemv(CblasNoTrans, 1.0, &subT, &subv, beta, &subp);
    /* last part of the product p */
    subp = gsl_vector_subvector(p, p->size - (i + 1) * D, D).vector;
    subv = gsl_vector_const_subvector(v, v->size - temp, temp).vector;
    subT = gsl_matrix_submatrix(T, 0, (s + i) * bCols - temp, D, temp).matrix;    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &subT, &subv, beta, &subp);
  }
  /* middle part */
  for (i = s_1; i < m - s_1; i++) {
    subp = gsl_vector_subvector(p, i * T->size1, T->size1).vector;
    subv = gsl_vector_const_subvector(v, (i - s_1) * bCols, T->size2).vector;
    gsl_blas_dgemv(CblasNoTrans, 1.0, T, &subv, beta, &subp);
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
