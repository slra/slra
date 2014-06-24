#include <memory.h>
#include <cstdarg>
#include "slra.h"

StationaryCholesky::StationaryCholesky( const StationaryStructure *s,  size_t d ) : 
                                      MuDependentCholesky(s, d), myStStruct(s)  {
  myGammaK = gsl_matrix_alloc(d, d * (getMu() + 1));
}  
  
StationaryCholesky::~StationaryCholesky() {
  gsl_matrix_free(myGammaK);
}

void StationaryCholesky::computeGammak( const gsl_matrix *Rt, double reg ) {
  gsl_matrix_view submat;
  
  for (size_t k = 0; k < getMu(); k++) { 
    submat = gsl_matrix_submatrix(myGammaK, 0, k * getD(), getD(), getD());
    myStStruct->AtVkB(&submat.matrix, k, Rt, Rt, myTempVijtRt);
 
    if (reg > 0) {
      gsl_vector diag = gsl_matrix_diagonal(&submat.matrix).vector;
      gsl_vector_add_constant(&diag, reg);
    }    
  }
  submat = gsl_matrix_submatrix(myGammaK, 0, getMu() * getD(), getD(), getD());
  gsl_matrix_set_zero(&submat.matrix);
}
  
void StationaryCholesky::computeGammaUpperTrg( const gsl_matrix *R, double reg ) {
  computeGammak(R, reg);
  
  size_t row_gam, col_gam, icor;
  double *gp = myPackedCholesky;
    
  for (size_t i = 0; i < myDMu; i++) {
    for (size_t j = 0; j < getD(); j++) {
      icor = i + j + 1;
      gp[i + j * myDMu] = gsl_matrix_get(myGammaK, 
          icor % getD(), j + (getMu() - (icor / getD())) * getD());
    }
  }
  for (size_t r = 1; r < getN(); r++) {
    gp +=  myDMu * getD();
    memcpy(gp, myPackedCholesky, myDMu * getD() * sizeof(double));
  }
}


