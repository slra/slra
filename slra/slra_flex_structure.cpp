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


slraGammaComputations *slraFlexStructure::createGammaComputations( int r, double reg_gamma )  {
  return new slraFlexGammaComputations(this, r, getM() / getK(), 1, reg_gamma);
}

slraDerivativeComputations *slraFlexStructure::createDerivativeComputations( int r ) {
  return new slraFlexDerivativeComputations(this, r, getK());
}

  


slraFlexStructure::slraFlexStructure( const slraFlexStructure &s ) :
    myK(s.myK), myQ(s.myQ), mySA(NULL) {
  int k;  
  
  mySA = new slraFlexBlock[myQ];
  for (k = 0; k < getQ(); k++) { 
    mySA[k] = s.mySA[k];
  }
  
  computeStats(); 
  computeWkParams();
  setNp(s.myNp);
}



slraFlexStructure::slraFlexStructure( const double *s_matr, size_t q, size_t k, int s_matr_cols, int np_or_m, bool set_m, 
                                      const double *w_k  ) :
                                      myQ(q), myK(k), mySA(NULL)  {
  mySA = new slraFlexBlock[myQ];
  for (size_t l = 0; l < myQ; l++) {
    mySA[l].blocks_in_row = *(s_matr + l);
    mySA[l].nb = (s_matr_cols > 1) ? *(s_matr + myQ + l): 1;
    
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

slraFlexStructure::~slraFlexStructure() {
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


slraFlexStructure::slraFlexStructure( const data_struct *s, int np ) : myQ(s->q), myK(s->k), mySA(NULL) {
  mySA = new slraFlexBlock[myQ];

  for (int l = 0; l < myQ; l++) {
    mySA[l]  = s->a[l];
  }
  
  computeStats();
  computeWkParams();
 
  setNp(np);                                      
}

void slraFlexStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  {
  int sum_np = 0, sum_nl = 0;
  int m_div_k = getM() / getK();
  size_t l,j, k, L;
  gsl_matrix_view c_chunk, c_chunk_sub;
 
  for (l = 0; l < getQ(); l++) {
    L = getFlexBlockLag(l);

    gsl_matrix_const_view p_matr_chunk = gsl_matrix_const_view_array(&(p->data[sum_np]), 
					 getFlexBlockT(l), getFlexBlockNb(l) * getK());

    for (k = 0; k < getK(); k++) {		   
      c_chunk = gsl_matrix_submatrix(c, k * m_div_k, sum_nl, m_div_k, getFlexBlockNCol(l));
				   
      for (j = 0; j < L; j++) {
        gsl_matrix_const_view p_matr_chunk_sub = gsl_matrix_const_submatrix(&p_matr_chunk.matrix, 
	                                        j, k * getFlexBlockNb(l), m_div_k, getFlexBlockNb(l));
        c_chunk_sub = gsl_matrix_submatrix(&c_chunk.matrix, 0,
                         j * getFlexBlockNb(l), m_div_k, getFlexBlockNb(l));
        gsl_matrix_memcpy(&c_chunk_sub.matrix, &p_matr_chunk_sub.matrix);
      }
      
    }  
    sum_np += getFlexBlockNp(l);
    sum_nl += getFlexBlockNCol(l);
  }
}

void slraFlexStructure::computeWkParams() {
  int k, l, i, offset, imax, sum_nl;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;
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
      offset = getFlexBlockNb(l) * k;
      imax   = ncol - offset;
      for (i = 0; i < imax; i++) {
        gsl_matrix_set(&zkl.matrix, i + offset, i, 1);
      }
      gsl_matrix_scale(&zkl.matrix, getInvBlockWeight(l));

      sum_nl += ncol;
    }
    myA[k] = zk;
  }
}

void slraFlexStructure::computeStats() {
  myNplusD = 0;
  myNpScale = 0;
  myNpOffset = 0;
  myMaxLag = 1;
  
  for (int l = 0; l < myQ; l++) {
    myNplusD += getFlexBlockNCol(l);
    myNpOffset += myK * (getFlexBlockLag(l) - 1) * getFlexBlockNb(l);
    myNpScale += getFlexBlockNb(l);
    
    if ((!isFlexBlockExact(l)) && getFlexBlockLag(l) > myMaxLag) {
      myMaxLag = mySA[l].blocks_in_row;
    }
  }
}


void slraFlexStructure::setNp( int np ) {
  if (np <= 0) {
    np = myNpScale + myNpOffset;
  }
  
  if (np < myNpScale + myNpOffset) {
    throw slraException("There is no matrix with the structure specification " 
                        "for length %d: vector too short.\n", np);
  }
  if ((np - myNpOffset) % myNpScale != 0) {
    throw slraException("There is no matrix with the structure specification "
                        "for length %d. scale = %d, offset = %d\n", np,
	                myNpScale, myNpOffset);
  }

  myNp = np;
}


void slraFlexStructure::correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
  int l, j, k, i, L, T;
  int sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix_view  p_matr_chunk, p_matr_chunk_sub;
  gsl_matrix_view brgf_matr, res_sub_matr, b_xext;
  gsl_vector_view p_chunk_vec, brgf_matr_row, res_sub;
  gsl_vector *res;
  double tmp;

  brgf_matr = gsl_matrix_view_vector(yr, getM(), R->size2);
    
  /* Allocate a vector for intermediate results */  
  res = gsl_vector_alloc(R->size1);

  for (l = 0; l < getQ(); l++) {
    int ncol = getFlexBlockNCol(l);
    b_xext = gsl_matrix_submatrix(R, sum_nl, 0, ncol, R->size2); 
    
    /* Calculate dimensions of current parameter chunk */
    L = getFlexBlockLag(l);
    T = getFlexBlockT(l);
    p_len = getFlexBlockNp(l);

    /* Adjust the result size */
    res_sub = gsl_vector_subvector(res, 0, ncol);
    res_sub_matr = gsl_matrix_view_vector(&res_sub.vector, L, getFlexBlockNb(l));
    /* Subtract correction if needed */
    if (!isFlexBlockExact(l)) {
      p_matr_chunk = gsl_matrix_view_array(&(p->data[sum_np]), T, getFlexBlockNb(l) * getK());
      for (j = 0; j < getK(); j++) {
        for (k = 0; k < (getM() / getK()); k++) {
          p_matr_chunk_sub = gsl_matrix_submatrix(&p_matr_chunk.matrix, 
				 k, getFlexBlockNb(l) * j, L, getFlexBlockNb(l));
          brgf_matr_row = gsl_matrix_row(&brgf_matr.matrix, k + j * (getM() / getK())); 
          gsl_blas_dgemv(CblasNoTrans, 1.0, &b_xext.matrix, 
			 &brgf_matr_row.vector, 0.0, &res_sub.vector);
	  gsl_vector_scale(&res_sub.vector, sqrt(getInvBlockWeight(l)));
          gsl_matrix_sub(&p_matr_chunk_sub.matrix, &res_sub_matr.matrix); 
        }
      }
    }
    sum_np += p_len;
    sum_nl += ncol;
  }
  
  gsl_vector_free(res);
}

void slraFlexStructure::setM( int m ) {
  if (m <= 0) {
    m = 1;
  }
  
  myNp = myNpScale * m + myNpOffset;
}



slraGammaComputations *slraFlexStructureExt::createGammaComputations( int r, double reg_gamma )  {
  return new slraFlexGammaComputationsExt(this, r, 1, reg_gamma);
}


slraDerivativeComputations *slraFlexStructureExt::createDerivativeComputations( int r ) {
  return new slraFlexDerivativeComputationsExt(this, r);
}



slraFlexStructureExt::slraFlexStructureExt( size_t q, size_t N, double *oldNk, double *oldMl, double *Wk ) :
    mySimpleStruct(oldNk, q, 1, 1, -1, false, Wk), myN(N), myOldMl(NULL) {
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

slraFlexStructureExt::~slraFlexStructureExt()  {
  if (myOldMl != NULL) {
    delete[] myOldMl;
  }
/*  if (myWk != NULL) {
    delete[] myWk;
  }*/
}

void slraFlexStructureExt::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  {
  int n_row = 0, sum_np = 0;
  gsl_matrix_view sub_c;

  
  for (int k = 0; k < getBlocksN(); sum_np += mySimpleStruct.getNp(), n_row += getMl(k), k++) {
    sub_c = gsl_matrix_submatrix(c, n_row, 0, getMl(k), c->size2);    
    mySimpleStruct.setM(getMl(k));
    gsl_vector_const_view sub_p = gsl_vector_const_subvector(p, sum_np, mySimpleStruct.getNp());
  
    mySimpleStruct.fillMatrixFromP(&sub_c.matrix, &sub_p.vector);
  }
}



void slraFlexStructureExt::correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) {
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




