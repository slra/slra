#include "slra.h"

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
}

/* Cholesky striped classes */
typedef Cholesky* pGammaCholesky;

StripedCholesky::StripedCholesky( const StripedStructure *s, size_t D ) : myS(s) {
  myD = D;
  myNGamma = myS->isSameGamma() ? 1 : myS->getBlocksN();
  myGamma = new pGammaCholesky[myNGamma];
  
  if (myNGamma == 1) {
    myGamma[0] = myS->getMaxBlock()->createCholesky(D); 
  } else {
    for (size_t k = 0; k < myS->getBlocksN(); k++) {
      myGamma[k] = myS->getBlock(k)->createCholesky(D);
    }
  }
}    

StripedCholesky::~StripedCholesky() {
  if (myGamma != NULL) {
    for (size_t k = 0; k < myNGamma; k++) {
      if (myGamma[k] != NULL) {
        delete myGamma[k];
      }
    }
    delete[] myGamma;
  }
}

void StripedCholesky::calcGammaCholesky( const gsl_matrix *R, double reg_gamma ) {
  for (size_t k = 0; k < myNGamma; k++) {
    myGamma[k]->calcGammaCholesky(R, reg_gamma);  
  }
}

void StripedCholesky::multInvCholeskyVector( gsl_vector * yr, long trans ) {
  size_t n_row = 0, k;
  gsl_vector yr_b;
  
  for (k = 0; k < myS->getBlocksN(); n_row += myS->getBlock(k)->getN(), k++) {
    yr_b = gsl_vector_subvector(yr, n_row*myD, 
                                myS->getBlock(k)->getN()*myD).vector;    
    myGamma[myNGamma == 1 ? 0 : k]->multInvCholeskyVector(&yr_b, trans);
  }
}

void StripedCholesky::multInvGammaVector( gsl_vector * yr ) {
  size_t n_row = 0, k;
  gsl_vector yr_b;
  
  for (k = 0; k < myS->getBlocksN(); n_row += myS->getBlock(k)->getN(), k++) {
    yr_b = gsl_vector_subvector(yr, n_row*myD, 
                                myS->getBlock(k)->getN()*myD).vector;    
    myGamma[myNGamma == 1 ? 0 : k]->multInvGammaVector(&yr_b);
  }
}

