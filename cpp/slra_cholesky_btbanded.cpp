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

StationaryCholesky::StationaryCholesky( const StationaryStructure *s,  int D ) : 
                                      SDependentCholesky(s, D), myWs(s)  {
  myGamma = gsl_matrix_alloc(getD(), getD() * (getS() + 1));
  myWkTmp = gsl_matrix_alloc(getM(), getD());
}  
  
StationaryCholesky::~StationaryCholesky() {
  gsl_matrix_free(myGamma);
  gsl_matrix_free(myWkTmp);
}

void StationaryCholesky::computeGammak( const gsl_matrix *R, double reg ) {
  gsl_matrix_view submat;
  
  for (size_t k = 0; k < getS(); k++) { /* compute brgamma_k = R' * w_k * R */
    submat = gsl_matrix_submatrix(myGamma, 0, k * getD(), getD(), getD());
    myWs->AtWkB(&submat.matrix, k, R, R, myWkTmp);
 
    if (reg > 0) {
      gsl_vector diag = gsl_matrix_diagonal(&submat.matrix).vector;
      gsl_vector_add_constant(&diag, reg);
    }    
  }
  submat = gsl_matrix_submatrix(myGamma, 0, getS() * getD(), getD(), getD());
  gsl_matrix_set_zero(&submat.matrix);
}
  
void StationaryCholesky::computeGammaUpperPart( const gsl_matrix *R, double reg ) {
  computeGammak(R, reg);
  
  int row_gam, col_gam, icor;
  double *gp = myPackedCholesky;
    
  for (int i = 0; i < d_times_s; i++) {
    for (int j = 0; j < getD(); j++) {
      icor = i + j + 1;
      gp[i + j * d_times_s] = gsl_matrix_get(myGamma, 
          icor % getD(), j + (getS() - (icor / getD())) * getD());
    }
  }
  for (int r = 1; r < getN(); r++) {
    gp +=  d_times_s * getD();
    memcpy(gp, myPackedCholesky, d_times_s * getD() * sizeof(double));
  }
}

SameStripedStationaryCholesky::
    SameStripedStationaryCholesky( const MosaicHStructure *s, 
         int D, int use_slicot  ) :  myS(s) {
  myBase = (StationaryCholesky *)myS->getMaxBlock()->createCholesky(D);  
}
SameStripedStationaryCholesky::~SameStripedStationaryCholesky() {
  delete myBase;
}
  
void SameStripedStationaryCholesky::calcGammaCholesky( const gsl_matrix *R, double reg_gamma ) {
  myBase->calcGammaCholesky(R, reg_gamma);
}

  
void SameStripedStationaryCholesky::
         multInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0, k;
  gsl_vector_view sub_yr;
  
  for (k = 0; k < myS->getBlocksN(); n_row += myS->getBlock(k)->getN(), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), 
                 myS->getBlock(k)->getN() * myBase->getD());    
  
    myBase->multInvPartCholeskyArray(sub_yr.vector.data, trans, 
        sub_yr.vector.size, sub_yr.vector.size);
  }
}

void SameStripedStationaryCholesky::multInvGammaVector( gsl_vector * yr ) {
  int n_row = 0, k;
  gsl_vector_view sub_yr;
  
  for (k = 0; k < myS->getBlocksN(); n_row += myS->getBlock(k)->getN(), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), 
                 myS->getBlock(k)->getN() * myBase->getD());    
    myBase->multInvPartGammaArray(sub_yr.vector.data, 
        sub_yr.vector.size, sub_yr.vector.size);
  }
}
