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

void StripedStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) {
  size_t n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;
  
  for (size_t k = 0; k < getBlocksN(); 
       sum_np += myStripe[k]->getNp(), n_row += getBlock(k)->getN(), k++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getBlock(k)->getN(), c->size2);    
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, 
        myStripe[k]->getNp());
    myStripe[k]->fillMatrixFromP(&sub_c.matrix, &sub_p.vector);
  }
}

void StripedStructure::multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
         const gsl_vector *y, double alpha, double beta, bool skipFixedBlocks ){
  size_t n_row = 0, sum_np = 0, D = R->size2;
  gsl_vector subp, suby;
  
  for (size_t k = 0; k < getBlocksN(); 
       sum_np += myStripe[k]->getNp(), n_row += getBlock(k)->getN() * D, k++) {
    suby = gsl_vector_const_subvector(y, n_row, getBlock(k)->getN() * D).vector;    
    subp = gsl_vector_subvector(p, sum_np, myStripe[k]->getNp()).vector;
    myStripe[k]->multByGtUnweighted(&subp, Rt, &suby, alpha, beta, 
                                    skipFixedBlocks);
  }                         
}

void StripedStructure::multByWInv( gsl_vector* p, long deg ) {
  size_t sum_np = 0;
  gsl_vector sub_p;
  
  for (size_t k = 0; k < getBlocksN(); sum_np += myStripe[k]->getNp(), k++) {
    sub_p = gsl_vector_subvector(p, sum_np, myStripe[k]->getNp()).vector;
    myStripe[k]->multByWInv(&sub_p, deg);
  }                         
}

Cholesky *StripedStructure::createCholesky( size_t D ) const {
  return new StripedCholesky(this, D);
}

DGamma *StripedStructure::createDGamma( size_t D ) const {
  return new StripedDGamma(this, D);
}




