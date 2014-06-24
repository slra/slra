#include "slra.h"

MuDependentCholesky::MuDependentCholesky( const MuDependentStructure *s,
                                          size_t d ) : myStruct(s), myD(d) {
  /* Calculate variables for FORTRAN routines */     
  myMu_1 =  myStruct->getMu() - 1;   // Maximal block superdiagonal
  myDMu =  myD * myStruct->getMu();  // 
  myDN = myStruct->getN() * myD;
  myDMu_1 = myD * myStruct->getMu() - 1;
  /* Preallocate arrays */
  myPackedCholesky = (double*)malloc(myDN * myDMu * sizeof(double));
  myTempVijtRt = gsl_matrix_alloc(myStruct->getM(), myD);
  myTempGammaij = gsl_matrix_alloc(myD, myD);
}
  
MuDependentCholesky::~MuDependentCholesky() {
  free(myPackedCholesky);
  gsl_matrix_free(myTempVijtRt);
  gsl_matrix_free(myTempGammaij);
}

void MuDependentCholesky::multInvCholeskyVector( gsl_vector * y_r, long trans ) {
  if (y_r->stride != 1) {
    throw new Exception("Cannot multiply vectors with stride != 1\n");
  }
  if (y_r->size > myDN) {
    throw new Exception("y_r->size > d * n\n");
  }
  size_t one = 1, info;
  dtbtrs_("U", (trans ? "T" : "N"), "N", &y_r->size, &myDMu_1, &one, 
	        myPackedCholesky, &myDMu, y_r->data, &y_r->size, &info);
}

void MuDependentCholesky::multInvGammaVector( gsl_vector * y_r ) {
  if (y_r->stride != 1) {
    throw new Exception("Cannot multiply vectors with stride != 1\n");
  }
  if (y_r->size > myDN) {
    throw new Exception("y_r->size > d * n\n");
  }
  size_t one = 1, info;
  dpbtrs_("U", &y_r->size, &myDMu_1, &one, 
          myPackedCholesky, &myDMu, y_r->data, &y_r->size, &info);  
}

void MuDependentCholesky::calcGammaCholesky( const gsl_matrix *Rt, double reg ) {
  size_t info = 0;
  computeGammaUpperTrg(Rt);
  dpbtrf_("U", &myDN, &myDMu_1, myPackedCholesky, &myDMu, &info);
  if (info && reg > 0) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, "Gamma is singular (DPBTRF info = %d), "
        "adding regularization, reg = %f.\n", info, reg);
    computeGammaUpperTrg(Rt, reg);
    dpbtrf_("U", &myDN, &myDMu_1, myPackedCholesky, &myDMu, &info);
  }
  if (info) {
    throw new Exception("Gamma is singular (DPBTRF info = %d).\n", info); 
  }
}

void MuDependentCholesky::computeGammaUpperTrg( const gsl_matrix *Rt, double reg ) {
  gsl_matrix gamma_ij;
  gsl_vector diag;
  double *blockColPtr =  myPackedCholesky;
  for (size_t i = 0; i < getN();      // Iterate by block columns of myPackedCholesky
                    ++i, blockColPtr += getMu() * getD() * getD()) {
    if (getMu() > 1) {  // We can process using GSL
      gsl_matrix blk_row =   // Create a view of i-th block row  (transposed, because in GSL row-major order)
          gsl_matrix_view_array_with_tda(blockColPtr + myDMu_1, // Start from the first element on the main diagonal
              (getMu() + 1) * getD(), getD(), myDMu_1).matrix;  // Adjust TDA to go from dpbtrf ordering to normal
      for (size_t k = 0; (k <= getMu()) && (k < getN() - i); ++k) { // k-th block diagonal 
        gamma_ij = gsl_matrix_submatrix(&blk_row, k * getD(),      // Select k-th block  
                                        0, getD(), getD()).matrix;
        if (k < getMu()) { 
 
          myStruct->AtVijB(myTempGammaij, i+k, i, // Block from lower block triangle of Gamma, 
                           Rt, Rt, myTempVijtRt);     // because in GSL the 
        } else {                                      // row-major order is used.
          gsl_matrix_set_zero(myTempGammaij);    // Put zeros in the mu-th block
        }
        if (k == 0) { // If we are on the main block diagonal
          if (reg > 0) { // Add regularization if needed
            diag = gsl_matrix_diagonal(myTempGammaij).vector;
            gsl_vector_add_constant(&diag, reg);
          }
          copyLowerTrg(&gamma_ij, myTempGammaij);      // Should copy only triangle
        } else {                                       // in order to avoid overwriting
          gsl_matrix_memcpy(&gamma_ij, myTempGammaij); // Otherwise copy as matrix
        }
      }
    } else { // If only one block diagonal
      myStruct->AtVijB(myTempGammaij, i, i, Rt, Rt, myTempVijtRt);  
      if (reg > 0) { // Add regularization if needed
        diag = gsl_matrix_diagonal(myTempGammaij).vector;
        gsl_vector_add_constant(&diag, reg);
      }
      gamma_ij = gsl_matrix_view_array(blockColPtr, getD(), getD()).matrix;
      shiftLowerTrg(&gamma_ij, myTempGammaij);
    }
  }
}
