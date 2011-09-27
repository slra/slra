/*********************************
 * Header file fo STLS package
 *********************************/

/* stls.h: STLS header file */
#ifndef _STLS_H_
#define _STLS_H_

#if defined(BUILD_R_PACKAGE)

#include <R.h>
#define PRINTF Rprintf

#elif defined(MEX_OCTAVE) ||  defined(MEX_MATLAB)

#include "mex.h"
#define PRINTF mexPrintf

#else

#include <stdio.h>
#define PRINTF printf

#endif



#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h> /* Levenberge-Marquardt */
#include <gsl/gsl_multimin.h>      /* BFGS Newton-type     */

/* size of the work array for mb02gd */
#define LDWORK 1000
#define EITER 1 /* maximum number of iterations reached */

/* optimization options and output information structure */
typedef struct {
  /* input options */
  int maxiter, disp; /* displayed information: 1 - notify, 2 - final, 3 - iter, 4 - off */
  double epsrel, epsabs, epsgrad;
  /* output information */
  int iter;
  double fmin;
  double time;
} opt_and_info;

/* structure in the data matrix C = [ A B ] */ 
#define MAXQ 10	/* maximum number of blocks in C */
typedef struct {
  int k;	/* = rowdim(block in T/H blocks) */ 
  int q;	/* number of blocks in C = [C1 ... Cq] */
  struct {
    char type;	/* 'T'-Toeplitz, 'H'-Hankel, 
                   'U'-unstructured, 'E'-exact */
    int ncol;	/* number of columns */
    int nb;     /* number of columns in a block of T/H blocks */
  } a[MAXQ];	/* q-element array describing C1,...,Cq; */
} data_struct;

/* three dimensional array of covariances w */
typedef struct {
  int s;	/* length of the array (w.s = s+1 from the paper) */
  gsl_matrix **a;
} w_data;

/* data needed for cost function and Jacobian evaluation */
typedef struct {
  gsl_matrix* a;
  gsl_matrix* b;
  w_data* w;
  int k; 
  int n_plus_d, 		/* = col_dim(C) */
    n_times_d,			/* = number of elements in x */
    k_times_d,			/* = row_dim(gamma) */
    k_times_d_times_s,		/* = col_dim(gamma) */
    k_times_d_times_s_minus_1,  /* = col_dim(gamma) - 1 */
    m_times_d, 			/* = row_dim(rb) */
    m_div_k, s_minus_1, 
    size_of_gamma, size_of_rb;
  /* Preallocated arrays */  
  gsl_matrix *x_ext; 
  double *rb;   /* Result of Cholesky factorization */
  gsl_vector *yr;
  
  /* Preallocated arrays for cholgam */
  gsl_matrix *tmp; /* Temp matrix for cholgam (x_ext' * w_k) P->k_times_d x SIZE_W  */
  gsl_matrix *gamma;
  double *gamma_vec;
  int ldwork;       /* Size of Dwork for MB02GD  */
  double *dwork;    /* Dwork for MB02GD  */
  
  /* Preallocated arrays for jacobian */
  gsl_matrix *dgamma, *st;
  double *jres1, * jres2;

  
  
} stls_opt_data;

/* Prototypes of functions */

#ifdef __cplusplus
extern "C" {
#endif

int stls(gsl_matrix*, gsl_matrix*, const data_struct*, 
	 gsl_matrix*, gsl_matrix*, opt_and_info* );
int stls_f (const gsl_vector*, void*, gsl_vector*);
double stls_f_ (const gsl_vector*, void*);
int stls_df (const gsl_vector*, void*, gsl_matrix*);
int stls_fdf (const gsl_vector*, 
	      void*, gsl_vector*, gsl_matrix*);
void print_state (int, gsl_multifit_fdfsolver*);
int s2w(const data_struct*, w_data*);
void print_mat(const gsl_matrix*);
void print_arr(double*, int);
void gsl_matrix_vectorize(double*, gsl_matrix*);
void gsl_matrix_vec_inv(gsl_matrix*, double*);
void tmv_prod(gsl_matrix*, int, 
	      gsl_vector*, int, gsl_vector*);
int tls(gsl_matrix*, gsl_matrix*, gsl_matrix*);

 
/* TODO: replace with something that uses printf. #define print_vec(v)   gsl_vector_fprintf(stdout,v,"%16.14f") */



void xmat2xext( gsl_matrix_const_view, gsl_matrix*, stls_opt_data* );
void cholgam( stls_opt_data* );
void jacobian( stls_opt_data*,  gsl_matrix*);

/* SLICOT and LAPACK functions */
/*
void mb02gd_(char*, char*, int*, int*, int*, const int*, int*, double*, int*, double*, int*, double*, const int*, int*);
void mb02md_(char*, int*, int*, int*, const int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, const int*, int*);
void dtbtrs_(char*, char*, char*, int*, int*, const int*, double*, int*, double*, int*, int*);
void dpbtrs_(char*, int*, int*, const int*, double*, int*, double*, int*, int*);
*/


void m_to_gsl_matrix(gsl_matrix* a_gsl, double* a_m);
void gsl_to_m_matrix(double* a_m, gsl_matrix* a_gsl); 


#ifdef __cplusplus
}
#endif


#endif /* _STLS_H_ */



