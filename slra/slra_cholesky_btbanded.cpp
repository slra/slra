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

StationaryCholesky::StationaryCholesky( const StationaryStructure *s, 
    int D, double reg_gamma  ) :  SDependentCholesky(s, D, reg_gamma), myWs(s)  {
  myGamma = gsl_matrix_alloc(getD(), getD() * (getS() + 1));
  myWkTmp = gsl_matrix_alloc(getNplusD(), getD());
}  
  
StationaryCholesky::~StationaryCholesky() {
  gsl_matrix_free(myGamma);
  gsl_matrix_free(myWkTmp);
}

void StationaryCholesky::computeGammak( gsl_matrix *R ) {
  size_t k;
  gsl_matrix_view submat;
  
  for (k = 0; k < getS(); k++) { /* compute brgamma_k = R' * w_k * R */
    submat = gsl_matrix_submatrix(myGamma, 0, k * getD(), getD(), getD());
    myWs->AtWkB(&submat.matrix, k, R, R, myWkTmp);
  }
  submat = gsl_matrix_submatrix(myGamma, 0, getS() * getD(), getD(), getD());
  gsl_matrix_set_zero(&submat.matrix);
}
  
void StationaryCholesky::computeGammaUpperPart( gsl_matrix *R ) {
  computeGammak(R);
  
  int row_gam, col_gam, icor;
  double *gp = myPackedCholesky;
    
  for (int i = 0; i < d_times_s; i++) {
    for (int j = 0; j < getD(); j++) {
      icor = i + j + 1;
      gp[i + j * d_times_s] = gsl_matrix_get(myGamma, 
          icor % getD(), j + (getS() - (icor / getD())) * getD());
    }
  }
  for (int r = 1; r < getM(); r++) {
    gp +=  d_times_s * getD();
    memcpy(gp, myPackedCholesky, d_times_s * getD() * sizeof(double));
  }
}

SameStripedStationaryCholesky::
    SameStripedStationaryCholesky( const MosaicHStructure *s, 
         int D, int use_slicot, double reg_gamma  ) :  myStruct(s) {
  myBase = (StationaryCholesky *)myStruct->getMaxBlock()->createCholesky(D, reg_gamma);  
}
SameStripedStationaryCholesky::~SameStripedStationaryCholesky() {
  delete myBase;
}
  
void SameStripedStationaryCholesky::calcGammaCholesky( gsl_matrix *R ) {
  myBase->calcGammaCholesky(R);
}

  
void SameStripedStationaryCholesky::
         multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), 
                 myStruct->getMl(k) * myBase->getD());    
  
    myBase->multiplyInvPartCholeskyArray(sub_yr.vector.data, trans, 
        sub_yr.vector.size, sub_yr.vector.size);
  }
}

void SameStripedStationaryCholesky::multiplyInvGammaVector( gsl_vector * yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myBase->getD(), 
                 myStruct->getMl(k) * myBase->getD());    
    myBase->multiplyInvPartGammaArray(sub_yr.vector.data, 
        sub_yr.vector.size, sub_yr.vector.size);
  }
}


