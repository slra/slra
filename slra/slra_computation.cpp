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

slraFlexComputationsParams::slraFlexComputationsParams( const slraFlexStructure *s ) {
  int k, l, i, offset, imax, sum_nl;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;
  char err_msg[70];
  int rep;
  int ncol;
  
  slraFlexComputationsParams *w = this;

  myS = s->getMaxLag();
  myA = (gsl_matrix**) malloc(myS * sizeof(gsl_matrix *));
  
  /* construct w */
  for (k = 0; k < myS; k++) { 
    zk   = gsl_matrix_alloc(s->getNplusD(), s->getNplusD());
    gsl_matrix_set_zero(zk);
    sum_nl = 0;
    for (l = 0; l < s->getQ(); l++) { 
      ncol = s->getFlexBlockNCol(l);
      zkl = gsl_matrix_submatrix(zk, sum_nl, sum_nl, ncol, ncol); 

      if ((!s->isFlexBlockExact(l))) {
	offset = s->getFlexBlockNb(l) * k;
	imax   = ncol - offset;
	for (i = 0; i < imax; i++) {
	  if (s->isFlexBlockToeplitz(l)) {
	    gsl_matrix_set(&zkl.matrix, i, i + offset, 1);
          } else {
	    gsl_matrix_set(&zkl.matrix, i + offset, i, 1);
	  }
	}
      } 
      
      sum_nl += ncol;
    }
    myA[k] = zk;
  }
}

slraFlexComputationsParams::~slraFlexComputationsParams() {
  for (int k = 0; k < myS; k++) {
    gsl_matrix_free(myA[k]);
  }
  free(myA);
}


slraFlexGammaComputations::slraFlexGammaComputations( int k, 
    int m, int n, int d, int use_slicot, slraFlexComputationsParams *w  ) : 
     myK(k), myM(m), myN(n), myD(d), my_use_slicot(use_slicot)  {
     
  myW = w;
    
  /* Preallocate arrays */
  myGammaVec = (double*) malloc(myD * myD * (myW->getS() + 1) * sizeof(double));
  myGamma = gsl_matrix_alloc(myD, myD * (myW->getS() + 1));
  myWkTmp = gsl_matrix_alloc(myD, myN + myD);
  myPackedCholesky = (double *)malloc((myM / myK) * myD * myD * myW->getS() * 
                              sizeof(double));
  myCholeskyWorkSize = 1 + myW->getS() * myD * myD + /* pDW */ 
                       3 * myD + /* 3 * K */
                       mymax(myW->getS(), (myM / myK) - myW->getS()) * myD * myD;
  myCholeskyWork = (double *)malloc(myCholeskyWorkSize * sizeof(double));                       
                       
                         
  /* Calculate variables for FORTRAN routines */     
  m_div_k = myM / myK;
  s_minus_1 = myW->getS() - 1;
  d_times_s = myD * myW->getS();
  d_times_m_div_k = (myM / myK) * myD;
  d_times_s_minus_1 = myD * myW->getS() - 1;
}
  
slraFlexGammaComputations::~slraFlexGammaComputations() {
  free(myGammaVec);
  gsl_matrix_free(myGamma);
  gsl_matrix_free(myWkTmp);
  free(myPackedCholesky);
  free(myCholeskyWork);
}
  
void slraFlexGammaComputations::computeCholeskyOfGamma( gsl_matrix *R )  {
  int k, info;
  gsl_matrix_view submat;
  const int zero = 0;

  /* compute brgamma_k = R' * w_k * R */
  for (k = 0; k < myW->getS(); k++) {
    submat = gsl_matrix_submatrix(myGamma, 0, k * myD, myD, myD);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, R, myW->getWk(k), 0.0, myWkTmp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myWkTmp, R, 0.0, 
                   &submat.matrix);
  }
  submat = gsl_matrix_submatrix(myGamma, 0, myW->getS() * myD, myD, myD);
  gsl_matrix_set_zero(&submat.matrix);
    
  if (my_use_slicot) { /* use SLICOT */
    gsl_matrix_vectorize(myGammaVec, myGamma);
    /* Cholesky factorization of Gamma */
    mb02gd_("R", "N", &myD, &m_div_k, &s_minus_1, &zero, 
            &m_div_k, myGammaVec, &myD, myPackedCholesky, &d_times_s, 
            myCholeskyWork, &myCholeskyWorkSize, &info); /**/
  } else { /* use LAPACK */
    int i, j, r, row_gam, col_gam, icor;
    double *gp = myPackedCholesky;
    
    for (i = 0; i < d_times_s; i++) {
      for (j = 0; j < myD; j++) {
        icor = i + j + 1;
        gp[i + j * d_times_s] = gsl_matrix_get(myGamma, 
            icor % myD, j + (myW->getS() - (icor / myD)) * myD);
      }
    }
    
    for (r = 1; r < m_div_k; r++) {
      gp +=  d_times_s * myD;
      memcpy(gp, myPackedCholesky, d_times_s * myD * sizeof(double));
    }
    
    dpbtrf_("U", &d_times_m_div_k, &d_times_s_minus_1, myPackedCholesky, 
            &d_times_s, &info);
  }
  
  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED */
  }
}
  
void slraFlexGammaComputations::multiplyInvCholeskyArray( double * yr, int trans, int rep ) {
  int info;
  int total_cols = myK * rep;

  dtbtrs_("U", (trans ? "T" : "N"), "N", 
          &d_times_m_div_k, &d_times_s_minus_1, &total_cols, 
	  myPackedCholesky, &d_times_s, yr, &d_times_m_div_k, &info);
}
  
void slraFlexGammaComputations::multiplyInvGammaArray( double * yr ) {
  int info;
  
  dpbtrs_("U", &d_times_m_div_k, &d_times_s_minus_1, &myK, 
          myPackedCholesky, &d_times_s, yr, &d_times_m_div_k, &info);  
}

void slraFlexGammaComputations::correctVector( gsl_vector* p, slraFlexStructure *s, 
                                               gsl_matrix *R, gsl_vector *f ) {
 
  int l, j, k, i, L, T;
  int sum_np = 0, sum_nl = 0, p_len;
  gsl_matrix_view  p_matr_chunk, p_matr_chunk_sub;
  gsl_matrix_view brgf_matr, res_sub_matr, b_xext;
  gsl_vector_view p_chunk_vec, brgf_matr_row, res_sub, b_xext_sub;
  gsl_matrix *xext_rev;
  gsl_vector *res;
  double tmp;

  brgf_matr = gsl_matrix_view_vector(f, myM, myD);
    
  /* Create reversed Xext matrix for Toeplitz blocks */
  xext_rev = gsl_matrix_alloc(myN + myD, myD); 

  for (i = 0; i < myN + myD; i++) {
    b_xext_sub = gsl_matrix_row(R, myN + myD -i-1);
    gsl_matrix_set_row(xext_rev, i,  &b_xext_sub.vector);
  }

  /* Allocate a vector for intermediate results */  
  res = gsl_vector_alloc(myN + myD);

  for (l = 0; l < s->getQ(); l++) {
    int ncol = s->getFlexBlockNCol(l);
  
    /* Select submatrix being used */
    if (s->isFlexBlockToeplitz(l)) {
      b_xext = gsl_matrix_submatrix(xext_rev, (myN + myD - (sum_nl + ncol)), 0, ncol, myD);
    } else {
      b_xext = gsl_matrix_submatrix(R, sum_nl, 0, ncol, myD); 
    }

    /* Calculate dimensions of current parameter chunk */
    L = s->getFlexBlockLag(l);
    T = s->getFlexBlockT(l);
    p_len = s->getFlexBlockNp(l);

    /* Adjust the result size */
    res_sub = gsl_vector_subvector(res, 0, ncol);
    res_sub_matr = gsl_matrix_view_vector(&res_sub.vector, L, s->getFlexBlockNb(l));
    /* Subtract correction if needed */
    if (!s->isFlexBlockExact(l)) {
      p_matr_chunk = gsl_matrix_view_array(&(p->data[sum_np]), T, s->getFlexBlockNb(l) * myK);
      for (j = 0; j < myK; j++) {
        for (k = 0; k < (myM / myK); k++) {
          p_matr_chunk_sub = gsl_matrix_submatrix(&p_matr_chunk.matrix, 
				 k, s->getFlexBlockNb(l) * j, L, s->getFlexBlockNb(l));
          brgf_matr_row = gsl_matrix_row(&brgf_matr.matrix, k + j * (myM / myK)); 
          gsl_blas_dgemv(CblasNoTrans, 1.0, &b_xext.matrix, 
			 &brgf_matr_row.vector, 0.0, &res_sub.vector);
          gsl_matrix_sub(&p_matr_chunk_sub.matrix, &res_sub_matr.matrix); 
        }
      }
    }
    sum_np += p_len;
    sum_nl += ncol;
  }
  
  gsl_vector_free(res);
  gsl_matrix_free(xext_rev); 
}



slraFlexDerivativeComputations::slraFlexDerivativeComputations( int k, 
    int m, int n, int d,  slraFlexComputationsParams *w  ) :
    myK(k), myM(m), myN(n), myD(d) {
  myW = w;
  d_times_m_div_k = (myM / myK) * myD;
  m_div_k = (myM / myK);
  
  myTempWkColRow = gsl_vector_alloc(myN + myD);
  myDGamma = gsl_matrix_alloc(myD, myD * (2 * myW->getS() - 1));
  
  mySubPhiT_Wk_R =  gsl_matrix_alloc(myN, myD);
  mySubPhiT_WkT_R = gsl_matrix_alloc(myN, myD);
  mySubPhiTmp  = gsl_matrix_alloc(myN + myD, myD);
  myN_k = gsl_matrix_alloc(myD, myD);
}

slraFlexDerivativeComputations::~slraFlexDerivativeComputations() {
  gsl_vector_free(myTempWkColRow);
  gsl_matrix_free(myDGamma);
  gsl_matrix_free(mySubPhiT_Wk_R);
  gsl_matrix_free(mySubPhiT_WkT_R);
  gsl_matrix_free(mySubPhiTmp);
  gsl_matrix_free(myN_k);
}


void slraFlexDerivativeComputations::computeYrtDgammaYr( gsl_matrix *mgrad, 
         gsl_matrix *R, gsl_matrix *perm, gsl_vector *yr ) {
  int k, ik;
  gsl_vector_view yr_sub, grad_sub;
  gsl_matrix_view yr_sub_matr, yr_sub_matr1, 
                  yr_sub_matr2, phi_sub_matr;
  phi_sub_matr = gsl_matrix_submatrix(perm, 0, 0, myN + myD, myN);

  gsl_matrix_set_zero(mgrad);
         
  for (k = 0; k < myW->getS(); k++) {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 2.0, myW->getWk(k), R, 0.0, 
                   mySubPhiTmp);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &phi_sub_matr.matrix, 
		   mySubPhiTmp, 0.0, mySubPhiT_Wk_R);

    if (k > 0) {
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, myW->getWk(k), R, 0.0, 
                     mySubPhiTmp);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &phi_sub_matr.matrix, 
		     mySubPhiTmp, 0.0, mySubPhiT_WkT_R);
    }


    for (ik = 0; ik < myK; ik++) {
      yr_sub = gsl_vector_subvector(yr, ik * d_times_m_div_k, d_times_m_div_k);
      yr_sub_matr = gsl_matrix_view_vector(&yr_sub.vector, m_div_k, myD);
      yr_sub_matr1 = gsl_matrix_submatrix(&yr_sub_matr.matrix, 
					  0, 0, m_div_k - k, myD);
      yr_sub_matr2 = gsl_matrix_submatrix(&yr_sub_matr.matrix, 
                                          k, 0, m_div_k - k, myD);
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, 
		     &yr_sub_matr1.matrix, &yr_sub_matr2.matrix, 
		     0.0, myN_k);

      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, mySubPhiT_Wk_R, myN_k, 1.0, mgrad);
      if (k > 0) {
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mySubPhiT_WkT_R, myN_k, 1.0, mgrad);
      }
    }
  }       
         
}

void slraFlexDerivativeComputations::computeDijGammaYr( gsl_vector *res, 
         gsl_matrix *R, gsl_matrix *perm, int i, int j, gsl_vector *yr ) {
  gsl_matrix_view dgamma_k, dgamma_minus_k;
  gsl_vector_view tmp1_row, tmp1_col;
  gsl_vector_view w_k_row, w_k_col;
  int k, l;
     
  /* form dgamma = d/d x_ij gamma */
  gsl_matrix_set_zero(myDGamma);

  for (k = 0; k < myW->getS(); k++) {
    dgamma_k = gsl_matrix_submatrix(myDGamma, 0, 
                   (myW->getS() - 1) * myD + k * myD, myD, myD);

    /* compute tmp1 = dx_ext' * w_k * x_ext * /
    /* Iterate over rows of dx_ext' * w_k */
    tmp1_row = gsl_matrix_row (&dgamma_k.matrix, j);
    w_k_row = gsl_matrix_column (perm, i);
    gsl_blas_dgemv(CblasTrans, 1.0, myW->getWk(k), &w_k_row.vector,
                   0.0, myTempWkColRow); 
    gsl_blas_dgemv(CblasTrans, 1.0, R, myTempWkColRow, 0.0, &tmp1_row.vector); 

    /* compute submat = submat  + x_ext' * tmp * dx_ext * /
    /* Iterate over rows of dx_ext' * w_k' */
    tmp1_col = gsl_matrix_column (&dgamma_k.matrix, j);
    gsl_blas_dgemv(CblasNoTrans, 1.0, myW->getWk(k), &w_k_row.vector,
                   0.0, myTempWkColRow); 
    gsl_blas_dgemv(CblasTrans, 1.0, R, myTempWkColRow, 1.0, &tmp1_col.vector); 
  }

  for (l = 0; l < myW->getS() - 1; l++) {
    dgamma_minus_k = gsl_matrix_submatrix(myDGamma, 0, l * myD, myD, myD);
    dgamma_k = gsl_matrix_submatrix(myDGamma, 0, 
                   (2 * myW->getS() - 2 - l) * myD, myD, myD);
    gsl_matrix_transpose_memcpy(&dgamma_minus_k.matrix, &dgamma_k.matrix);
  }
       
  /* compute st_ij = DGamma * yr */
  for (l = 0; l < myK; l++) { 
    gsl_vector_view res_v = gsl_vector_subvector(res, 
                                l * d_times_m_div_k, d_times_m_div_k);
    gsl_vector_view yr_v = gsl_vector_subvector(yr, 
                               l * d_times_m_div_k, d_times_m_div_k);
    tmv_prod_new(myDGamma, myW->getS(), &yr_v.vector, m_div_k, &res_v.vector);  
  }
}


slraException::slraException( const char *format, ... ) { 
  va_list vl;
  va_start(vl, format);  
  myMsg[MSG_MAX-1] = 0;
  vsnprintf(myMsg, MSG_MAX-1, format, vl); 
}


slraFlexStructure::slraFlexStructure( const double *s_matr, int q, int k, int s_matr_cols, int np ) :
                                      myQ(q), myK(k)  {
  for (int l = 0; l < myQ; l++) {
    mySA[l].blocks_in_row = *(s_matr + l);
    mySA[l].nb = (s_matr_cols > 1) ? *(s_matr + myQ + l): 1;
    mySA[l].exact = (s_matr_cols > 2) ? *(s_matr + 2 * myQ + l): 0;
    mySA[l].toeplitz = (s_matr_cols > 3) ? *(s_matr + 3 * myQ + l): 0;
  }    
   
  computeStats(); 
  setNp(np);                                      
}

slraFlexStructure::slraFlexStructure( const data_struct *s, int np ) : myQ(s->q), myK(s->k) {
  for (int l = 0; l < myQ; l++) {
    mySA[l]  = s->a[l];
  }
  
  computeStats(); 
  setNp(np);                                      
}

void slraFlexStructure::fillMatrixFromP( gsl_matrix* c, gsl_vector* p ) const {
  int sum_np = 0, sum_nl = 0;
  int m_div_k = getM() / getK();
  int l,j, k, L;
  gsl_matrix_view p_matr_chunk, c_chunk, p_matr_chunk_sub, c_chunk_sub;
 
  for (l = 0; l < getQ(); l++) {
    L = getFlexBlockLag(l);

    p_matr_chunk = gsl_matrix_view_array(&(p->data[sum_np]), 
					 getFlexBlockT(l), getFlexBlockNb(l) * getK());

    for (k = 0; k < getK(); k++) {		   
      c_chunk = gsl_matrix_submatrix(c, k * m_div_k, sum_nl, m_div_k, getFlexBlockNCol(l));
				   
      for (j = 0; j < L; j++) {
        p_matr_chunk_sub = gsl_matrix_submatrix(&p_matr_chunk.matrix, 
	                                        j, k * getFlexBlockNb(l), m_div_k, getFlexBlockNb(l));
        c_chunk_sub = gsl_matrix_submatrix(&c_chunk.matrix, 0, 
                      (isFlexBlockToeplitz(l) ? (L- j -1) * getFlexBlockNb(l) : j * getFlexBlockNb(l)), 
                      m_div_k, getFlexBlockNb(l));
        gsl_matrix_memcpy(&c_chunk_sub.matrix, &p_matr_chunk_sub.matrix);
      }
      
    }  
    sum_np += getFlexBlockNp(l);
    sum_nl += getFlexBlockNCol(l);
  }
}



void slraFlexStructure::computeStats() {
  myNplusD = 0;
  myNpScale = 0;
  myNpOffset = 0;
  myMaxLag = 1;
  
  for (int l = 0; l < myQ; l++) {
    myNplusD += mySA[l].blocks_in_row * mySA[l].nb;
    myNpOffset += myK * (mySA[l].blocks_in_row - 1) * mySA[l].nb;
    myNpScale += mySA[l].nb;
    
    if ((!mySA[l].exact) && mySA[l].blocks_in_row > myMaxLag) {
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





/*


class slraFnComputation {


  void computeLsPseudoJacobianFromYr( gsl_matrix *lsJacobian, const gsl_matrix *R, 
  {
    computePseudoJacobianFirstTerm(lsJacobian);
    computePseudoJacobianSecondTerm(P->tmpJacobian2, R, Phi);
    
    gsl_matrix_scale(P->tmpJacobian2, 0.5);
    gsl_matrix_sub(lsJacobian, P->tmpJacobian2);
  }


  void computeLsFun( gsl_vector *lsFun, const gsl_matrix *R ) {
    computeCholeskyOfGamma(R);
    computeYr(LsFunValue, R);
    multiplyInvCholesky(LsFunValue, 1);
  }


  void computeLsPseudoJacobianOnly( gsl_matrix *lsJacobian, const gsl_matrix *R, 
           const gsl_matrix *Phi ) {
    computeCholeskyOfGamma(R);
    computeYr(P->tmpYr, R);

    computeLsPseudoJacobianFromYr(computeLsPseudoJacobianFromYr, R, Phi);
  }

  void computeLsPseudoJacobianAndFun( gsl_vector *lsFun, gsl_matrix *lsJacobian, 
           const gsl_matrix *R, const gsl_matrix *Phi ) {
    computeCholeskyOfGamma(R);
    computeYr(P->tmpYr, R);
    
       
  }
}*/

