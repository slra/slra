#include "slra.h"

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
}

/* Structure striped classes */
StripedStructure::StripedStructure( size_t N, slraStructure **stripe  ): myN(N), myStripe(stripe) {
  size_t k;  
  
  for (k = 0, myM = 0, myNp = 0, myMaxMlInd = 0; k < myN; 
       myM += myStripe[k]->getM(), myNp += myStripe[k]->getNp(), k++) {

    if (myStripe[k]->getM() > myStripe[myMaxMlInd]->getM()) {
      myMaxMlInd = k;
    }
  }
}

StripedStructure::~StripedStructure()  {
  if (myStripe != NULL) {
    for (size_t k = 0; k < myN; k++) {
      if (myStripe[k] != NULL) {
        delete myStripe[k];
      }
    }
    delete[] myStripe;
  }
}

void StripedStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  {
  int n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;
  
  for (int k = 0; k < getBlocksN(); sum_np += myStripe[k]->getNp(), n_row += getMl(k), k++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getMl(k), c->size2);    
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, myStripe[k]->getNp());
    myStripe[k]->fillMatrixFromP(&sub_c.matrix, &sub_p.vector);
  }
}

void StripedStructure::correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  size_t n_row = 0, sum_np = 0, D = R->size2;
  gsl_vector_view sub_p, sub_yr;
  
  
  for (int k = 0; k < getBlocksN(); 
       sum_np += myStripe[k]->getNp(), n_row += getMl(k) * D, k++) {
    PRINTF("Computing correction:\n");

    sub_yr = gsl_vector_subvector(yr, n_row, getMl(k) * D);    
    sub_p = gsl_vector_subvector(p, sum_np, myStripe[k]->getNp());
    DEBUGINT(myStripe[k]);
    myStripe[k]->correctVector(&sub_p.vector, R, &sub_yr.vector);
  }
}

slraGammaCholesky *StripedStructure::createGammaComputations( int r, double reg_gamma ) const {
  return new StripedCholesky(this, r, reg_gamma);
}

slraDGamma *StripedStructure::createDerivativeComputations( int r ) const {
  return new StripedDGamma(this, r);
}

/* Cholesky striped classes */
typedef slraGammaCholesky* pslraGammaCholesky;

StripedCholesky::StripedCholesky( const StripedStructure *s, int r, double reg_gamma  ) :
    myStruct(s) {
  myD = s->getNplusD() - r;
  myGamma = new pslraGammaCholesky[myStruct->getBlocksN()];
  for (int k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k] = myStruct->getBlock(k)->createGammaComputations(r, reg_gamma);
  }
}    

StripedCholesky::~StripedCholesky() {
  if (myGamma != NULL) {
    for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
      if (myGamma[k] != NULL) {
        delete myGamma[k];
      }
    }
    delete[] myGamma;
  }
}

void StripedCholesky::multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
    myGamma[k]->multiplyInvCholeskyVector(&sub_yr.vector, trans);
  }
}

void StripedCholesky::multiplyInvGammaVector( gsl_vector * yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
    myGamma[k]->multiplyInvGammaVector(&sub_yr.vector);
  }
}

void StripedCholesky::computeCholeskyOfGamma( gsl_matrix *R ) {
  for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k]->computeCholeskyOfGamma(R);  
  }
}

/* DGamma striped classes */
typedef slraDGamma* pslraDGamma;

StripedDGamma::StripedDGamma( const StripedStructure *s, int r  ) : 
    myStruct(s) {
  myTmpGrad = gsl_matrix_alloc(myStruct->getNplusD(), myStruct->getNplusD() - r);  
  myLHDGamma = new pslraDGamma[myStruct->getBlocksN()];
  for (int k = 0; k < myStruct->getBlocksN(); k++) {
    myLHDGamma[k] = myStruct->getBlock(k)->createDerivativeComputations(r);
  }
}    

StripedDGamma::~StripedDGamma() {
  gsl_matrix_free(myTmpGrad);
  if (myLHDGamma != NULL) {
    for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
      if (myLHDGamma[k] != NULL) {
        delete myLHDGamma[k];
      }
    }
    delete[] myLHDGamma;
  }
}

void StripedDGamma::calcYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, 
                                     gsl_vector *yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * R->size2, myStruct->getMl(k) * R->size2);    
    myLHDGamma[k]->calcYrtDgammaYr(myTmpGrad, R, &sub_yr.vector);
    gsl_matrix_add(grad, myTmpGrad);
  }
}

void StripedDGamma::calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                        gsl_matrix *perm, int i, int j, gsl_vector *yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr, sub_res;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * R->size2, myStruct->getMl(k) * R->size2);    
    sub_res = gsl_vector_subvector(res, n_row * R->size2, myStruct->getMl(k) * R->size2);    
    myLHDGamma[k]->calcDijGammaYr(&sub_res.vector, R, perm, i, j, &sub_yr.vector);
  }                   
}
