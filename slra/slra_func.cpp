#include <memory.h>
#include <time.h>

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

/*
 * tmv_prod: block-Toeplitz matrix T times vector v
 *
 * tt - storage for [t_s-1' ... t_1' t_0 t_1 ... t_s-1].
 * s - number of blocks in t, t = [t_0 ... t_s-1]
 * s_1 = s - 1;  m = (int) v->size1 / tt->size1
 * p - result
 */ 
void tmv_prod_new( gsl_matrix *tt, int s, gsl_vector* v, int m, 
		   gsl_vector* p)
{
  int i, imax, temp, s_1 = s - 1;
  int row_lim = GSL_MIN(s_1, m/2);
  gsl_vector_view subv;
  gsl_vector_view subp; 	/* subvector of v and p */

  int TM = tt->size1; 		/* = block size */

  gsl_matrix_view submat, source;

  /* construct p = T*v */
  gsl_vector_set_zero(p);

  /* beginning and end parts of the product p */
  for (i = 0; i < row_lim; i++) {
    temp = GSL_MIN(s+i, m)*TM;
    /* beginning part */
    subp = gsl_vector_subvector(p, i*TM, TM);
    subv = gsl_vector_subvector(v, 0, temp);
    submat = gsl_matrix_submatrix
      (tt, 0, (s_1-i)*TM, TM, temp);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &submat.matrix, 
		   &subv.vector, 0.0, &subp.vector);
    /* last part */
    subp = gsl_vector_subvector(p, p->size - (i+1)*TM, TM);
    subv = gsl_vector_subvector(v, v->size - temp, temp);
    submat = gsl_matrix_submatrix(tt, 0, (s+i)*TM -temp, TM, temp);    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &submat.matrix, 
		   &subv.vector, 0.0, &subp.vector);
  }

  /* middle part */
  for (i = s_1, imax = m - s_1 ; i < imax; i++) {
    subp = gsl_vector_subvector(p, i*TM, TM);
    subv = gsl_vector_subvector(v, (i-s_1)*TM, tt->size2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, tt, &subv.vector, 
		   0.0, &subp.vector);
  }
  
}


