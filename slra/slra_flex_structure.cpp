#include <limits>


#include <memory.h>
#include <cstdarg>

extern "C" {

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

}

#include "slra.h"

slraException::slraException( const char *format, ... ) { 
  va_list vl;
  va_start(vl, format);  
  myMsg[MSG_MAX-1] = 0;
  vsnprintf(myMsg, MSG_MAX-1, format, vl); 
}


slraGammaComputations *slraLayeredHankelStructure::createGammaComputations( int r, double reg_gamma )  {
  return new slraFlexGammaComputations(this, r, getM(), 1, reg_gamma);
}

slraDerivativeComputations *slraLayeredHankelStructure::createDerivativeComputations( int r ) {
  return new slraFlexDerivativeComputations(this, r);
}

  


slraLayeredHankelStructure::slraLayeredHankelStructure( const slraLayeredHankelStructure &s ) :  myQ(s.myQ), mySA(NULL) {
  int k;  
  
  mySA = new slraFlexBlock[myQ];
  for (k = 0; k < getQ(); k++) { 
    mySA[k] = s.mySA[k];
  }
  
  computeStats(); 
  computeWkParams();
  setNp(s.myNp);
}



slraLayeredHankelStructure::slraLayeredHankelStructure( const double *oldNk, size_t q, int np_or_m, bool set_m, const double *w_k  ) :
                                      myQ(q), mySA(NULL)  {
  mySA = new slraFlexBlock[myQ];
  for (size_t l = 0; l < myQ; l++) {
    mySA[l].blocks_in_row = oldNk[l];
    
    if (w_k != NULL) {
      if (w_k[l] != std::numeric_limits<double>::infinity() && w_k[l] <= 0) {
        throw new slraException("This value of weight is not supported: %lf\n", w_k[l]);
      }
      mySA[l].inv_w = 1 / w_k[l];
    } else {
      mySA[l].inv_w = 1;
    }
  }    
   
  computeStats(); 
  computeWkParams();
  if (set_m) {
    setM(np_or_m);                                      
  } else {
    setNp(np_or_m);                                      
  }
}

slraLayeredHankelStructure::~slraLayeredHankelStructure() {
  if (mySA != NULL) {
    delete[] mySA;
  }
  if (myA != NULL) {
    for (int k = 0; k < myMaxLag; k++) {
      gsl_matrix_free(myA[k]);
    }
    free(myA);
  }
}

void slraLayeredHankelStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  {
  int sum_np = 0, sum_nl = 0;
  size_t l, j;
  gsl_matrix_view c_chunk, c_chunk_sub;
 
  for (l = 0; l < getQ(); l++) {
    gsl_vector_const_view p_chunk = gsl_vector_const_subvector(p, sum_np, getFlexBlockNp(l));
    c_chunk = gsl_matrix_submatrix(c, 0, sum_nl, getM(), getFlexBlockNCol(l));
				   
    for (j = 0; j < getFlexBlockLag(l); j++) {
      gsl_vector_const_view p_chunk_sub = 
          gsl_vector_const_subvector(&p_chunk.vector, j, getM());
      gsl_matrix_set_col(&c_chunk.matrix, j, &p_chunk_sub.vector);
    }  
    sum_np += getFlexBlockNp(l);
    sum_nl += getFlexBlockNCol(l);
  }
}

void slraLayeredHankelStructure::computeWkParams() {
  int k, l, i, imax, sum_nl;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;
  gsl_vector_view diag;
  int rep, ncol;
  int myS = getMaxLag();

  myA = (gsl_matrix**) malloc(myS * sizeof(gsl_matrix *));
  
  /* construct w */
  for (k = 0; k < myS; k++) { 
    zk   = gsl_matrix_alloc(getNplusD(), getNplusD());
    gsl_matrix_set_zero(zk);
    sum_nl = 0;

    for (l = 0; l < getQ(); l++) { 
      ncol = getFlexBlockNCol(l);
      zkl = gsl_matrix_submatrix(zk, sum_nl, sum_nl, ncol, ncol); 
      if (k < ncol) {
        diag = gsl_matrix_subdiagonal(&zkl.matrix, k);
        gsl_vector_set_all(&diag.vector, getInvBlockWeight(l)); 
      }
      sum_nl += ncol;
    }
    myA[k] = zk;
  }
}

void slraLayeredHankelStructure::computeStats() {
  myNplusD = 0;
  myNpOffset = 0;
  myMaxLag = 1;
  
  for (int l = 0; l < myQ; l++) {
    myNplusD += getFlexBlockNCol(l);
    myNpOffset += getFlexBlockLag(l) - 1;
    
    if ((!isFlexBlockExact(l)) && getFlexBlockLag(l) > myMaxLag) {
      myMaxLag = mySA[l].blocks_in_row;
    }
  }
}


void slraLayeredHankelStructure::setNp( int np ) {
  if (np <= 0) {
    np = getNpScale() + myNpOffset;
  }
  
  if (np < getNpScale() + myNpOffset) {
    throw slraException("There is no matrix with the structure specification " 
                        "for length %d: vector too short.\n", np);
  }
  if ((np - myNpOffset) % getNpScale() != 0) {
    throw slraException("There is no matrix with the structure specification "
                        "for length %d. scale = %d, offset = %d\n", np,
	                getNpScale(), myNpOffset);
  }

  myNp = np;
}


void slraLayeredHankelStructure::correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  int l, k;
  int sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix_view brgf_matr, b_xext;
  gsl_vector_view brgf_matr_row, res_sub, p_chunk_sub;
  gsl_vector *res = gsl_vector_alloc(R->size1);

  brgf_matr = gsl_matrix_view_vector(yr, getM(), R->size2);

  for (l = 0; l < getQ(); l++) {
    b_xext = gsl_matrix_submatrix(R, sum_nl, 0, getFlexBlockNCol(l), R->size2); 
    res_sub = gsl_vector_subvector(res, 0, getFlexBlockNCol(l));
    if (!isFlexBlockExact(l)) {    /* Subtract correction if needed */
      for (k = 0; k < getM(); k++) {
        p_chunk_sub =  gsl_vector_subvector(p, k + sum_np, getFlexBlockLag(l));
        brgf_matr_row = gsl_matrix_row(&brgf_matr.matrix, k); 
        gsl_blas_dgemv(CblasNoTrans, sqrt(getInvBlockWeight(l)), &b_xext.matrix, 
                        &brgf_matr_row.vector, 0.0, &res_sub.vector); 
        gsl_vector_sub(&p_chunk_sub.vector, &res_sub.vector); 
      }
    }
    sum_np += getFlexBlockNp(l);
    sum_nl += getFlexBlockNCol(l);
  }
  
  gsl_vector_free(res);
}

void slraLayeredHankelStructure::setM( int m ) {
  if (m <= 0) {
    m = 1;
  }
  
  myNp = getNpScale() * m + myNpOffset;
}



slraGammaComputations *slraMosaicHankelStructure::createGammaComputations( int r, double reg_gamma )  {
  return new slraFlexGammaComputationsExt(this, r, 1, reg_gamma);
}


slraDerivativeComputations *slraMosaicHankelStructure::createDerivativeComputations( int r ) {
  return new slraFlexDerivativeComputationsExt(this, r);
}



slraMosaicHankelStructure::slraMosaicHankelStructure( size_t q, size_t N, double *oldNk, double *oldMl, double *Wk ) :
    mySimpleStruct(oldNk, q, -1, false, Wk), myN(N), myOldMl(NULL) {
  size_t k;  

  myOldMl = new size_t[N];  
  /*myWk = new double[q];
  
  for (k = 0; k < q; k++) {
    myWk[k] = (Wk != NULL ? Wk[k] : 1);
  }*/
  


  myM = 0;
  myMaxMl = 0;
  for (k = 0; k < myN; k++) {
    myOldMl[k] = oldMl[k];
    myM += myOldMl[k];
    
    if (myOldMl[k] > myMaxMl) {
      myMaxMl = myOldMl[k];
    }
  }
  
}

slraMosaicHankelStructure::~slraMosaicHankelStructure()  {
  if (myOldMl != NULL) {
    delete[] myOldMl;
  }
/*  if (myWk != NULL) {
    delete[] myWk;
  }*/
}

void slraMosaicHankelStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  {
  int n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;

  
  for (int k = 0; k < getBlocksN(); sum_np += mySimpleStruct.getNp(), n_row += getMl(k), k++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getMl(k), c->size2);    
    mySimpleStruct.setM(getMl(k));
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, mySimpleStruct.getNp());
  
    mySimpleStruct.fillMatrixFromP(&sub_c.matrix, &sub_p.vector);
  }
}



void slraMosaicHankelStructure::correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  int n_row = 0, sum_np = 0;
  gsl_vector_view sub_p;
  gsl_vector_view sub_yr;
  int D = R->size2;
  
  for (int k = 0; k < getBlocksN(); sum_np += mySimpleStruct.getNp(), n_row += getMl(k) * D, k++) {
    sub_yr = gsl_vector_subvector(yr, n_row, getMl(k) * D);    
    mySimpleStruct.setM(getMl(k));
    sub_p = gsl_vector_subvector(p, sum_np, mySimpleStruct.getNp());
  
    mySimpleStruct.correctVector(&sub_p.vector, R, &sub_yr.vector);
  }
}




