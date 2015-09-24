#include "slra.h"

/* DGamma striped classes */
typedef DGamma* pDGamma;

StripedDGamma::StripedDGamma( const StripedStructure *s, size_t d  ) : 
    myS(s) {
  myTmpGrad = gsl_matrix_alloc(myS->getM(), d);  
  myLHDGamma = new pDGamma[myS->getBlocksN()];
  for (size_t k = 0; k < myS->getBlocksN(); k++) {
    myLHDGamma[k] = myS->getBlock(k)->createDGamma(d);
  }
}    

StripedDGamma::~StripedDGamma() {
  gsl_matrix_free(myTmpGrad);
  if (myLHDGamma != NULL) {
    for (size_t k = 0; k < myS->getBlocksN(); k++) {
      if (myLHDGamma[k] != NULL) {
        delete myLHDGamma[k];
      }
    }
    delete[] myLHDGamma;
  }
}

void StripedDGamma::calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt, 
                                   const gsl_matrix *Yt ) {
  gsl_matrix_set_zero(At);
  for (size_t n_row = 0, k = 0; k < myS->getBlocksN(); 
                                n_row += myS->getBlock(k)->getN(), ++k) {
    gsl_matrix subYt = gsl_matrix_const_submatrix(Yt, n_row, 0, 
                           myS->getBlock(k)->getN(), Rt->size2).matrix;    
    myLHDGamma[k]->calcYtDgammaY(myTmpGrad, Rt, &subYt);
    gsl_matrix_add(At, myTmpGrad);
  }
}

void StripedDGamma::calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt, 
                        size_t j_1, size_t i_1, const gsl_vector *y,
                        const gsl_matrix *Phi ) {
  size_t n_row = 0, k;
  gsl_vector sub_y, sub_z;
  
  for (k = 0; k < myS->getBlocksN(); n_row += myS->getBlock(k)->getN(), ++k) {
    sub_y = gsl_vector_const_subvector(y, n_row * Rt->size2, 
                                  myS->getBlock(k)->getN() * Rt->size2).vector;    
    sub_z = gsl_vector_const_subvector(z, n_row * Rt->size2, 
                                   myS->getBlock(k)->getN() * Rt->size2).vector;    
    myLHDGamma[k]->calcDijGammaYr(&sub_z, Rt, j_1, i_1, &sub_y, Phi);
  }                   
}
