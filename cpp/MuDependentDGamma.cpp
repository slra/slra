#include <memory.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

MuDependentDGamma::MuDependentDGamma( const MuDependentStructure *s, size_t d ) :
   myD(d), myStruct(s) {
  myTmp1 = gsl_vector_alloc(myD);  
  myTmp2 = gsl_vector_alloc(myStruct->getM());  
  myTmp3 = gsl_vector_alloc(myStruct->getM());  
  myEye = gsl_matrix_alloc(myStruct->getM(), myStruct->getM());
  gsl_matrix_set_identity(myEye);
}

MuDependentDGamma::~MuDependentDGamma(){
  gsl_vector_free(myTmp1);
  gsl_vector_free(myTmp2);
  gsl_vector_free(myTmp3);
  gsl_matrix_free(myEye);
}

 void MuDependentDGamma::calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt, 
                   const gsl_matrix *Yt ) {
  size_t n = Yt->size1; long Mu = myStruct->getMu();
  gsl_vector y_i, y_j;
  gsl_matrix tmp2_v = gsl_matrix_view_vector(myTmp2, myStruct->getM(), 1).matrix,
             tmp3_v = gsl_matrix_view_vector(myTmp3, myStruct->getM(), 1).matrix;

  gsl_matrix_set_zero(At);
  for (size_t i = 0; i < n; ++i) {
    y_i = gsl_matrix_const_row(Yt, i).vector;
    
    size_t j = (i + 1 > Mu ? i - Mu + 1 : 0);
    for (; j < mymin(i + Mu, n); ++j) {
      y_j = gsl_matrix_const_row(Yt, j).vector;

      gsl_blas_dgemv(CblasNoTrans, 1.0, Rt, &y_i, 0.0, myTmp2);
      myStruct->VijB(&tmp3_v, j, i, &tmp2_v);
      gsl_blas_dger(2.0, myTmp3, &y_j, At);
    }
  }
}

void MuDependentDGamma::calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt, 
                   size_t j_1, size_t i_1, const gsl_vector *y )  {
  gsl_vector e_j = gsl_matrix_column(myEye, j_1).vector, yr_sub, z_sub;
  size_t k, l, n = y->size / myD; long Mu = myStruct->getMu();
  double tmp;

  gsl_vector_set_zero(z); 
  for (k = 0; k < n; k++)  {
    z_sub = gsl_vector_subvector(z, k * myD, myD).vector;
 
    l = (k + 1 > Mu ? k - Mu + 1 : 0);
    for (;  l < mymin(k + Mu, n); l++) {
      yr_sub = gsl_vector_const_subvector(y, l * myD, myD).vector;
      myStruct->AtVijV(myTmp1, k, l, Rt, &e_j, myTmp2);
      gsl_blas_daxpy(gsl_vector_get(&yr_sub, i_1), myTmp1, &z_sub);

      myStruct->AtVijV(myTmp1, l, k, Rt, &e_j, myTmp2);
      gsl_blas_ddot(myTmp1, &yr_sub, &tmp);
      gsl_vector_set(&z_sub, i_1, tmp + gsl_vector_get(&z_sub, i_1));
    }
  }
}
