#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

MuDependentDGamma::MuDependentDGamma( const MuDependentStructure *s, size_t D ) :
   myD(D), myW(s) {
  myTmp1 = gsl_vector_alloc(myD);  
  myTmp2 = gsl_vector_alloc(myW->getM());  
  myTmp3 = gsl_vector_alloc(myW->getM());  
  myEye = gsl_matrix_alloc(myW->getM(), myW->getM());
  gsl_matrix_set_identity(myEye);
}

MuDependentDGamma::~MuDependentDGamma(){
  gsl_vector_free(myTmp1);
  gsl_vector_free(myTmp2);
  gsl_vector_free(myTmp3);
  gsl_matrix_free(myEye);
}

 void MuDependentDGamma::calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                   const gsl_vector *yr ) {
  size_t n = yr->size / myD; long Mu = myW->getMu();
  gsl_matrix Y_r = gsl_matrix_const_view_vector(yr, n, myD).matrix;
  gsl_vector y_i, y_j;
  gsl_matrix tmp2_v = gsl_matrix_view_vector(myTmp2, myW->getM(), 1).matrix,
             tmp3_v = gsl_matrix_view_vector(myTmp3, myW->getM(), 1).matrix;

  gsl_matrix_set_zero(grad);
  for (size_t i = 0; i < n; i++) {
    y_i = gsl_matrix_row(&Y_r, i).vector;
    
    size_t j = (i + 1 > Mu ? i - Mu + 1 : 0);
    for (; j < mymin(i + Mu, n); j++) {
      y_j = gsl_matrix_row(&Y_r, j).vector;

      gsl_blas_dgemv(CblasNoTrans, 1.0, R, &y_i, 0.0, myTmp2);
      myW->VijB(&tmp3_v, j, i, &tmp2_v);
      gsl_blas_dger(2.0, myTmp3, &y_j, grad);
    }
  }
}

void MuDependentDGamma::calcDijGammaYr( gsl_vector *res, const gsl_matrix *R, 
                   size_t i, size_t j, const gsl_vector *Yr ) {
  gsl_vector perm_col = gsl_matrix_column(myEye, i).vector, yr_sub, res_sub;
  size_t k, l, n = Yr->size / myD; long Mu = myW->getMu();
  double tmp;

  gsl_vector_set_zero(res); 
  for (k = 0; k < n; k++)  {
    res_sub = gsl_vector_subvector(res, k * myD, myD).vector;
 
    l = (k + 1 > Mu ? k - Mu + 1 : 0);
    for (;  l < mymin(k + Mu, n); l++) {
      yr_sub = gsl_vector_const_subvector(Yr, l * myD, myD).vector;
      myW->AtVijV(myTmp1, k, l, R, &perm_col, myTmp2);
      gsl_blas_daxpy(gsl_vector_get(&yr_sub, j), myTmp1, &res_sub);

      myW->AtVijV(myTmp1, l, k, R, &perm_col, myTmp2);
      gsl_blas_ddot(myTmp1, &yr_sub, &tmp);
      gsl_vector_set(&res_sub, j, tmp + gsl_vector_get(&res_sub, j));
    }
  }
}
