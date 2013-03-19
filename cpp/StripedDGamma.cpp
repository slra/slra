#include "slra.h"

/* DGamma striped classes */
typedef DGamma* pDGamma;

StripedDGamma::StripedDGamma( const StripedStructure *s, size_t D  ) : 
    myS(s) {
  myTmpGrad = gsl_matrix_alloc(myS->getM(), D);  
  myLHDGamma = new pDGamma[myS->getBlocksN()];
  for (size_t k = 0; k < myS->getBlocksN(); k++) {
    myLHDGamma[k] = myS->getBlock(k)->createDGamma(D);
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

void StripedDGamma::calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                                     const gsl_vector *yr ) {
  size_t n_row = 0, k;
  
  gsl_matrix_set_zero(grad);
  for (k = 0; k < myS->getBlocksN(); n_row += myS->getBlock(k)->getN(), k++) {
    gsl_vector_const_view sub_yr = gsl_vector_const_subvector(yr, n_row * R->size2, 
                                  myS->getBlock(k)->getN() * R->size2);    
    myLHDGamma[k]->calcYrtDgammaYr(myTmpGrad, R, &sub_yr.vector);
    gsl_matrix_add(grad, myTmpGrad);
  }
}

void StripedDGamma::calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                        size_t i, size_t j, gsl_vector *yr ) {
  size_t n_row = 0, k;
  gsl_vector_view sub_yr, sub_res;
  
  for (k = 0; k < myS->getBlocksN(); n_row += myS->getBlock(k)->getN(), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * R->size2, 
                                  myS->getBlock(k)->getN() * R->size2);    
    sub_res = gsl_vector_subvector(res, n_row * R->size2, 
                                   myS->getBlock(k)->getN() * R->size2);    
    myLHDGamma[k]->calcDijGammaYr(&sub_res.vector, R, i, j, &sub_yr.vector);
  }                   
}
