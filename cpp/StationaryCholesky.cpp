#include <memory.h>
#include <cstdarg>
#include "slra.h"

StationaryCholesky::StationaryCholesky( const StationaryStructure *s,  size_t D ) : 
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
  
  size_t row_gam, col_gam, icor;
  double *gp = myPackedCholesky;
    
  for (size_t i = 0; i < myDS; i++) {
    for (size_t j = 0; j < getD(); j++) {
      icor = i + j + 1;
      gp[i + j * myDS] = gsl_matrix_get(myGamma, 
          icor % getD(), j + (getS() - (icor / getD())) * getD());
    }
  }
  for (size_t r = 1; r < getN(); r++) {
    gp +=  myDS * getD();
    memcpy(gp, myPackedCholesky, myDS * getD() * sizeof(double));
  }
}


