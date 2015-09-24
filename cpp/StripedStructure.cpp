#include "slra.h"

/* Structure striped classes */
StripedStructure::StripedStructure( size_t blocksN, Structure **stripe, 
                                    bool isSameGamma  ) :
    myBlocksN(blocksN), myStripe(stripe), myIsSameGamma(isSameGamma) {
  size_t l;  
  
  for (l = 0, myN = 0, myNp = 0, myMaxNlInd = 0; l < myBlocksN; 
       myN += myStripe[l]->getN(), myNp += myStripe[l]->getNp(), l++) {

    if (myStripe[l]->getN() > myStripe[myMaxNlInd]->getN()) {
      myMaxNlInd = l;
    }
  }
}

StripedStructure::~StripedStructure()  {
  if (myStripe != NULL) {
    for (size_t l = 0; l < myBlocksN; l++) {
      if (myStripe[l] != NULL) {
        delete myStripe[l];
      }
    }
    delete[] myStripe;
  }
}

void StripedStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) {
  size_t n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;
  
  for (size_t l = 0; l < getBlocksN(); 
       sum_np += myStripe[l]->getNp(), n_row += getBlock(l)->getN(), l++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getBlock(l)->getN(), c->size2);    
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, 
        myStripe[l]->getNp());
    myStripe[l]->fillMatrixFromP(&sub_c.matrix, &sub_p.vector);
  }
}

void StripedStructure::multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
         const gsl_vector *y, double alpha, double beta, bool skipFixedBlocks ){
  size_t n_row = 0, sum_np = 0, d = Rt->size2;
  gsl_vector subp, suby;
  
  for (size_t l = 0; l < getBlocksN(); 
       sum_np += myStripe[l]->getNp(), n_row += getBlock(l)->getN() * d, l++) {
    suby = gsl_vector_const_subvector(y, n_row, getBlock(l)->getN() * d).vector;    
    subp = gsl_vector_subvector(p, sum_np, myStripe[l]->getNp()).vector;
    myStripe[l]->multByGtUnweighted(&subp, Rt, &suby, alpha, beta, 
                                    skipFixedBlocks);
  }                         
}

void StripedStructure::multByWInv( gsl_vector* p, long deg ) const {
  size_t sum_np = 0;
  gsl_vector sub_p;
  
  for (size_t l = 0; l < getBlocksN(); sum_np += myStripe[l]->getNp(), l++) {
    sub_p = gsl_vector_subvector(p, sum_np, myStripe[l]->getNp()).vector;
    myStripe[l]->multByWInv(&sub_p, deg);
  }                         
}

Cholesky *StripedStructure::createCholesky( size_t d ) const {
  return new StripedCholesky(this, d);
}

DGamma *StripedStructure::createDGamma( size_t d ) const {
  return new StripedDGamma(this, d);
}




