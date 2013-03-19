#include "slra.h"

/* Structure striped classes */
StripedStructure::StripedStructure( size_t blocksN, Structure **stripe, 
                                    bool isSameGamma  ) :
    myBlocksN(blocksN), myStripe(stripe), myIsSameGamma(isSameGamma) {
  size_t k;  
  
  for (k = 0, myN = 0, myNp = 0, myMaxNkInd = 0; k < myBlocksN; 
       myN += myStripe[k]->getN(), myNp += myStripe[k]->getNp(), k++) {

    if (myStripe[k]->getN() > myStripe[myMaxNkInd]->getN()) {
      myMaxNkInd = k;
    }
  }
}

StripedStructure::~StripedStructure()  {
  if (myStripe != NULL) {
    for (size_t k = 0; k < myBlocksN; k++) {
      if (myStripe[k] != NULL) {
        delete myStripe[k];
      }
    }
    delete[] myStripe;
  }
}

void StripedStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p,
                                bool premultInvW ) {
  size_t n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;
  
  for (size_t k = 0; k < getBlocksN(); 
       sum_np += myStripe[k]->getNp(), n_row += getBlock(k)->getN(), k++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getBlock(k)->getN(), c->size2);    
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, 
        myStripe[k]->getNp());
    myStripe[k]->fillMatrixFromP(&sub_c.matrix, &sub_p.vector, premultInvW);
  }
}

void StripedStructure::correctP( gsl_vector* p, const gsl_matrix *R, 
                                 const gsl_vector *yr, long wdeg ) {
  size_t n_row = 0, sum_np = 0, D = R->size2;
  gsl_vector_view sub_p;
  
  for (size_t k = 0; k < getBlocksN(); 
       sum_np += myStripe[k]->getNp(), n_row += getBlock(k)->getN()*D, k++) {
    gsl_vector_const_view sub_yr = 
        gsl_vector_const_subvector(yr, n_row, getBlock(k)->getN() * D);    
    sub_p = gsl_vector_subvector(p, sum_np, myStripe[k]->getNp());
    myStripe[k]->correctP(&sub_p.vector, R, &sub_yr.vector, wdeg);
  }
}

Cholesky *StripedStructure::createCholesky( size_t D ) const {
  return new StripedCholesky(this, D);
}

DGamma *StripedStructure::createDGamma( size_t D ) const {
  return new StripedDGamma(this, D);
}




