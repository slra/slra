#include "slra.h"


extern "C" {

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

}

/* Structure striped classes */

slraGammaCholesky *slraStripedStructure::createGammaComputations( int r, double reg_gamma ) const {
  return new slraDiagGammaCholesky(this, r, reg_gamma);
}

slraDGamma *slraStripedStructure::createDerivativeComputations( int r ) const {
  return new slraDGammaStriped(this, r);
}


slraStripedStructure::slraStripedStructure( size_t N, slraStructure **stripe  ) : myN(N), myLHStripe(stripe) {
  size_t k;  
  myM = 0;
  myNp = 0;
  myMaxMlInd = 0;

  for (k = 0; k < myN; k++) {
    myM += myLHStripe[k]->getM();
    myNp += myLHStripe[k]->getNp();
    if (myLHStripe[k]->getM() > myLHStripe[myMaxMlInd]->getM()) {
      myMaxMlInd = k;
    }
  }
}


slraStripedStructure::~slraStripedStructure()  {
  if (myLHStripe != NULL) {
    for (size_t k = 0; k < myN; k++) {
      if (myLHStripe[k] != NULL) {
        delete myLHStripe[k];
      }
    }
    delete[] myLHStripe;
  }
}

void slraStripedStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  {
  int n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;
  
  for (int k = 0; k < getBlocksN(); sum_np += myLHStripe[k]->getNp(), n_row += getMl(k), k++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getMl(k), c->size2);    
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, myLHStripe[k]->getNp());
    myLHStripe[k]->fillMatrixFromP(&sub_c.matrix, &sub_p.vector);
  }
}

void slraStripedStructure::correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  int n_row = 0, sum_np = 0;
  gsl_vector_view sub_p, sub_yr;
  int D = R->size2;
  
  for (int k = 0; k < getBlocksN(); sum_np += myLHStripe[k]->getNp(), n_row += getMl(k) * D, k++) {
    sub_yr = gsl_vector_subvector(yr, n_row, getMl(k) * D);    
    sub_p = gsl_vector_subvector(p, sum_np, myLHStripe[k]->getNp());
    myLHStripe[k]->correctVector(&sub_p.vector, R, &sub_yr.vector);
  }
}


/* Cholesky striped classes */

typedef slraGammaCholesky* pslraGammaCholesky;

slraDiagGammaCholesky::slraDiagGammaCholesky( const slraStripedStructure *s, int r, double reg_gamma  ) :
    myStruct(s) {
  myD = s->getNplusD() - r;
  myGamma = new pslraGammaCholesky[myStruct->getBlocksN()];
  for (int k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k] = myStruct->getBlock(k)->createGammaComputations(r, reg_gamma);
  }
}    

slraDiagGammaCholesky::~slraDiagGammaCholesky() {
  if (myGamma != NULL) {
    for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
      if (myGamma[k] != NULL) {
        delete myGamma[k];
      }
    }
    delete[] myGamma;
  }
}

void slraDiagGammaCholesky::multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
    myGamma[k]->multiplyInvCholeskyVector(&sub_yr.vector, trans);
  }
}

void slraDiagGammaCholesky::multiplyInvGammaVector( gsl_vector * yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * myD, myStruct->getMl(k) * myD);    
    myGamma[k]->multiplyInvGammaVector(&sub_yr.vector);
  }
}

void slraDiagGammaCholesky::computeCholeskyOfGamma( gsl_matrix *R ) {
  for (size_t k = 0; k < myStruct->getBlocksN(); k++) {
    myGamma[k]->computeCholeskyOfGamma(R);  
  }
}

/* DGamma striped classes */

typedef slraDGamma* pslraDGamma;

slraDGammaStriped::slraDGammaStriped( const slraStripedStructure *s, int r  ) :
    myStruct(s) {
  myTmpGrad = gsl_matrix_alloc(myStruct->getNplusD(), myStruct->getNplusD() - r);  
  
  myLHDGamma = new pslraDGamma[myStruct->getBlocksN()];
  for (int k = 0; k < myStruct->getBlocksN(); k++) {
    myLHDGamma[k] = myStruct->getBlock(k)->createDerivativeComputations(r);
  }
}    

slraDGammaStriped::~slraDGammaStriped() {
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


void slraDGammaStriped::computeYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, gsl_vector *yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * R->size2, myStruct->getMl(k) * R->size2);    
    myLHDGamma[k]->computeYrtDgammaYr(myTmpGrad, R, &sub_yr.vector);
    gsl_matrix_add(grad, myTmpGrad);
  }
}


void slraDGammaStriped::computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *yr ) {
  int n_row = 0;
  gsl_vector_view sub_yr;
  gsl_vector_view sub_res;
  
  for (int k = 0; k < myStruct->getBlocksN(); n_row += myStruct->getMl(k), k++) {
    sub_yr = gsl_vector_subvector(yr, n_row * R->size2, myStruct->getMl(k) * R->size2);    
    sub_res = gsl_vector_subvector(res, n_row * R->size2, myStruct->getMl(k) * R->size2);    
    myLHDGamma[k]->computeDijGammaYr(&sub_res.vector, R, perm, i, j, &sub_yr.vector);
  }                   
}
