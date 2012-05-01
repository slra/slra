#include "slra.h"

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
}

/* Structure striped classes */
StripedStructure::StripedStructure( size_t N, Structure **stripe  ) :
    myN(N), myStripe(stripe) {
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

void StripedStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) {
  int n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;
  
  for (int k = 0; k < getBlocksN(); 
       sum_np += myStripe[k]->getNp(), n_row += getMl(k), k++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getMl(k), c->size2);    
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, 
        myStripe[k]->getNp());
    myStripe[k]->fillMatrixFromP(&sub_c.matrix, &sub_p.vector);
  }
}

void StripedStructure::correctP( gsl_vector* p, gsl_matrix *R, 
                                 gsl_vector *yr, bool scaled ) {
  size_t n_row = 0, sum_np = 0, D = R->size2;
  gsl_vector_view sub_p, sub_yr;
  
  for (int k = 0; k < getBlocksN(); 
       sum_np += myStripe[k]->getNp(), n_row += getMl(k) * D, k++) {
    sub_yr = gsl_vector_subvector(yr, n_row, getMl(k) * D);    
    sub_p = gsl_vector_subvector(p, sum_np, myStripe[k]->getNp());
    myStripe[k]->correctP(&sub_p.vector, R, &sub_yr.vector, scaled);
  }
}

Cholesky *StripedStructure::createCholesky( int D, double reg_gamma ) const {
  return new StripedCholesky(this, D, reg_gamma);
}

DGamma *StripedStructure::createDGamma( int D ) const {
  return new StripedDGamma(this, D);
}

/* Cholesky striped classes */
typedef Cholesky* pGammaCholesky;

StripedCholesky::StripedCholesky( const StripedStructure *s, int D, 
                                   double reg_gamma ) :
    myStruct(s) {
  myD = D;
  myGamma = new pGammaCholesky[myStruct->getBlocksN()];
  for (int k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k] = myStruct->getBlock(k)->createCholesky(D, reg_gamma);
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

void StripedCholesky::multInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0, k;
  gsl_vector_view sub_yr;
  
  for (k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
    myGamma[k]->multInvCholeskyVector(&sub_yr.vector, trans);
  }
}

void StripedCholesky::multInvGammaVector( gsl_vector * yr ) {
  int n_row = 0, k;
  gsl_vector_view sub_yr;
  
  for (k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
    myGamma[k]->multInvGammaVector(&sub_yr.vector);
  }
}

void StripedCholesky::calcGammaCholesky( gsl_matrix *R ) {
  for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k]->calcGammaCholesky(R);  
  }
}

/* DGamma striped classes */
typedef DGamma* pDGamma;

StripedDGamma::StripedDGamma( const StripedStructure *s, int D  ) : 
    myStruct(s) {
  myTmpGrad = gsl_matrix_alloc(myStruct->getNplusD(), D);  
  myLHDGamma = new pDGamma[myStruct->getBlocksN()];
  for (int k = 0; k < myStruct->getBlocksN(); k++) {
    myLHDGamma[k] = myStruct->getBlock(k)->createDGamma(D);
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
  int n_row = 0, k;
  gsl_vector_view sub_yr;
  
  for (k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * R->size2, 
                                  myStruct->getMl(k) * R->size2);    
    myLHDGamma[k]->calcYrtDgammaYr(myTmpGrad, R, &sub_yr.vector);
    gsl_matrix_add(grad, myTmpGrad);
  }
}

void StripedDGamma::calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                        gsl_matrix *perm, int i, int j, gsl_vector *yr ) {
  int n_row = 0, k;
  gsl_vector_view sub_yr, sub_res;
  
  for (k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * R->size2, 
                                  myStruct->getMl(k) * R->size2);    
    sub_res = gsl_vector_subvector(res, n_row * R->size2, 
                                   myStruct->getMl(k) * R->size2);    
    myLHDGamma[k]->calcDijGammaYr(&sub_res.vector, R, perm, i, j,
                                  &sub_yr.vector);
  }                   
}
