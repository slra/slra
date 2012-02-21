#include <memory.h>

extern "C" {

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

}


#include "slra.h"

slraFlexGammaComputations::slraFlexGammaComputations( const data_struct *s, int m, int n, int d, int use_slicot  ) : 
                                         myS(*s), myM(m), myN(n), myD(d)  {
  s2w(s, &myW, myN + myD, 0);  
    
  /* Preallocate arrays */
  myGammaVec = (double*) malloc(myD * myD * (myW.s + 1) * sizeof(double));
  myGamma = gsl_matrix_alloc(myD, myD * (myW.s + 1));
  myWkTmp = gsl_matrix_alloc(myD, myN + myD);
  myPackedCholesky = (double *)malloc((myM / myS.k) * myD * myD * myW.s * 
                              sizeof(double));
  myCholeskyWorkSize = 1 + myW.s * myD * myD + /* pDW */ 
                       3 * myD + /* 3 * K */
                       mymax(myW.s, (myM / myS.k) - myW.s) * myD * myD;
  myCholeskyWork = (double *)malloc(myCholeskyWorkSize * sizeof(double));                       
                       
                         
  /* Calculate variables for FORTRAN routines */     
  m_div_k = myM / myS.k;
  s_minus_1 = myW.s - 1;
  d_times_s = myD * myW.s;
  d_times_m_div_k = (myM / myS.k) * myD;
  d_times_s_minus_1 = myD * myW.s - 1;
}
  
slraFlexGammaComputations::~slraFlexGammaComputations() {
  for (int k = 0; k < myW.s; k++) {
    gsl_matrix_free(myW.a[k]);
  }
  free(myW.a);
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
  for (k = 0; k < myW.s; k++) {
    submat = gsl_matrix_submatrix(myGamma, 0, k * myD, myD, myD);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, R, myW.a[k], 0.0, myWkTmp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, myWkTmp, R, 0.0, 
                   &submat.matrix);
  }
  submat = gsl_matrix_submatrix(myGamma, 0, myW.s * myD, myD, myD);
  gsl_matrix_set_zero(&submat.matrix);
    
  if (use_slicot) { /* use SLICOT */
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
            icor % myD, j + (myW.s - (icor / myD)) * myD);
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
  
void slraFlexGammaComputations::multiplyInvCholesky ( double * yr, int trans, int rep ) {
  int info;
  int total_cols = myS.k * rep;

  dtbtrs_("U", (trans ? "T" : "N"), "N", 
          &d_times_m_div_k, &d_times_s_minus_1, &total_cols, 
	  myPackedCholesky, &d_times_s, yr, &d_times_m_div_k, &info);
}
  
void slraFlexGammaComputations::multiplyInvGamma ( double * yr ) {
  int info;
  
  dpbtrs_("U", &d_times_m_div_k, &d_times_s_minus_1, &myS.k, 
          myPackedCholesky, &d_times_s, yr, &d_times_m_div_k, &info);  
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

