#include "slra.h"

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
}

/* Cholesky striped classes */
typedef Cholesky* pGammaCholesky;

StripedCholesky::StripedCholesky( const StripedStructure *s, size_t d ) : myStruct(s), 
                     myD(d), myNGamma(s->isSameGamma() ? 1 : s->getBlocksN())  {
  myGamma = new pGammaCholesky[myNGamma];
  
  if (myNGamma == 1) {
    myGamma[0] = myStruct->getMaxBlock()->createCholesky(d); 
  } else {
    for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
      myGamma[k] = myStruct->getBlock(k)->createCholesky(d);
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

void StripedCholesky::calcGammaCholesky( const gsl_matrix *Rt, double reg ) {
  for (size_t k = 0; k < myNGamma; k++) {
    myGamma[k]->calcGammaCholesky(Rt, reg);  
  }
}

void StripedCholesky::multInvCholeskyVector( gsl_vector * y_r, long trans ) {
  size_t n_row = 0, k;
  gsl_vector yr_b;
  
  for (k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getBlock(k)->getN(), k++) {
    yr_b = gsl_vector_subvector(y_r, n_row * myD, 
                                myStruct->getBlock(k)->getN() * myD).vector;    
    myGamma[myNGamma == 1 ? 0 : k]->multInvCholeskyVector(&yr_b, trans);
  }
}

void StripedCholesky::multInvGammaVector( gsl_vector * y_r ) {
  size_t n_row = 0, k;
  gsl_vector yr_b;
  
  for (k = 0; k < myStruct->getBlocksN();
              n_row += myStruct->getBlock(k)->getN(), k++) {
    yr_b = gsl_vector_subvector(y_r, n_row * myD, 
                                myStruct->getBlock(k)->getN() * myD).vector;    
    myGamma[myNGamma == 1 ? 0 : k]->multInvGammaVector(&yr_b);
  }
}

