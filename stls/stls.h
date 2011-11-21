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
#define EITER 1 /* maximum number of iterations reached */

#ifdef __cplusplus
extern "C" {
#endif


/* optimization options and output information structure */
typedef struct {
  /* input options */
  int maxiter, disp; /* displayed information: 1 - notify, 2 - final, 3 - iter, 4 - off */
  double epsrel, epsabs, epsgrad;
  
  double reggamma; /* To be worked out */
  /* output information */
  int iter;
  double fmin;
  double time;
} opt_and_info;

/* structure in the data matrix C = [ A B ] */ 
#define MAXQ 10	/* maximum number of blocks in C */
typedef struct {
  /* Parameters given by the user */
  int m;        /* number of rows */       
        
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



#define mymax(a,b) ((a) > (b) ? (a) : (b)) 
#define mymin(a,b) ((a) < (b) ? (a) : (b))


#define COMMON_PARAMS  \
    gsl_matrix* a; \
    gsl_matrix* b; \
    double reggamma; \
    w_data w; \
    int k;  \
    int  n_plus_d, 		/* = col_dim(C) */ \
    n_times_d,			/* = number of elements in x */ \
    k_times_d,			/* = row_dim(gamma) */ \
    k_times_d_times_s,		/* = col_dim(gamma) */ \
    k_times_d_times_s_minus_1,  /* = col_dim(gamma) - 1 */ \
    m_times_d, 			/* = row_dim(rb) */ \
    m_div_k,  \
    s_minus_1;  \
    int one; /* One for blas routines */ 


#define PREPARE_COMMON_PARAMS(A, B, S, OPT, PP, isblock) \
  do {\
   int m = A->size1, n = A->size2, d = B->size2;		\
  /* set other parameters */\
  PP->a = A;\
  PP->b = B;\
  PP->reggamma = OPT->reggamma; \
  /* find Wk  */ \
  s2w(S, &PP->w, n+d, isblock);\
  PP->k = S->k;  \
  PP->n_plus_d = n + d;   \
  PP->n_times_d = n * d;   \
  PP->k_times_d = S->k * d; \
  PP->k_times_d_times_s = PP->k_times_d * PP->w.s;\
  PP->k_times_d_times_s_minus_1 = PP->k_times_d_times_s - 1;\
  PP->m_times_d = m * d;\
  PP->m_div_k = (int) m / S->k;\
  PP->s_minus_1 = PP->w.s - 1;\
  PP->one = 1;\
  \
  } while (0)





typedef struct {
  COMMON_PARAMS;

  /* Preallocated arrays */  
  gsl_matrix *x_ext; 
  gsl_vector *yr;
  
  /* Preallocated arrays for cholgam */
  double *rb;   /* Result of Cholesky factorization */
  gsl_matrix *tmp; /* Temp matrix for cholgam (x_ext' * w_k) P->k_times_d x SIZE_W  */
  gsl_matrix *gamma;
  double *gamma_vec;
  int ldwork;       /* Size of Dwork for MB02GD  */
  double *dwork;    /* Dwork for MB02GD  */

  /* Preallocated arrays for jacobian */
  gsl_matrix *dgamma, *st;
  double *jres1, * jres2;
} stls_opt_data_old;



/* data needed for cost function and Jacobian evaluation */
typedef struct {
  COMMON_PARAMS;

  int m, n, d;
  /* Preallocated arrays for cholgam (new) */
  int  d_times_s;		/* = col_dim(new gamma) */
  int  d_times_s_minus_1;  /* = col_dim(new gamma) - 1 */
  int  d_times_m_div_k;  /* */

  double *rb2;   /* Result of Cholesky factorization (test) */

  /* Preallocated data for block computations */
  gsl_matrix *bx_ext;  

  double *brg_rb;   /* Result of Cholesky factorization */
  gsl_matrix *brg_tmp; /* Temp matrix for cholgam (x_ext' * w_k) P->k_times_d x SIZE_W  */
  gsl_matrix *brg_gamma;
  gsl_matrix *brg_dgamma;
  gsl_matrix *brg_tdgamma;
  gsl_matrix *brg_st;
  double *brg_jres2;
  
  double *brg_gamma_vec;
  int brg_ldwork;       /* Size of Dwork for MB02GD  */
  double *brg_dwork;    /* Dwork for MB02GD  */

  gsl_vector *brg_f;    /* Reshaped vector holding f */
  gsl_vector *brg_yr;    /* Reshaped vector holding y_r */
  
  gsl_matrix *brg_a, *brg_b;
  
  double *brg_j1b_vec;
  gsl_matrix *brg_j1b;
} stls_opt_data_reshaped;



/* Prototypes of functions */


int stls(gsl_matrix*, gsl_matrix*, data_struct*, 
	 gsl_matrix*, gsl_matrix*, opt_and_info*, gsl_vector *, int, int );
	 

int check_and_adjust_parameters( data_struct *s, int *n_plus_d, int *np );
int stls_fill_matrix_from_p( gsl_matrix* c,  data_struct *s, gsl_vector* p);
int stls_correction_reshaped(gsl_vector* p, data_struct *s, void* params, const gsl_vector* x);




void print_state (int, gsl_multifit_fdfsolver*);

int get_bandwidth_from_structure(const data_struct*);

int s2w(const data_struct*, w_data*, int,  int);
void print_mat(const gsl_matrix*);
void print_mat_tr(const gsl_matrix*);
void print_arr(double*, int);
void gsl_matrix_vectorize(double*, gsl_matrix*);
void gsl_matrix_vec_inv(gsl_matrix*, double*);
void tmv_prod(gsl_matrix*, int, 
	      gsl_vector*, int, gsl_vector*);
void tmv_prod_new(gsl_matrix*, int, 
	      gsl_vector*, int, gsl_vector*);

int tls(gsl_matrix*, gsl_matrix*, gsl_matrix*);

 
/* TODO: replace with something that uses printf. #define print_vec(v)   gsl_vector_fprintf(stdout,v,"%16.14f") */



void xmat2_block_of_xext( gsl_matrix_const_view, gsl_matrix *);


void allocate_and_prepare_data_reshaped( gsl_matrix* a, gsl_matrix* b, const data_struct* s, opt_and_info *opt, stls_opt_data_reshaped *P );
void free_memory_reshaped( stls_opt_data_reshaped *P );

double stls_f_reshaped_ (const gsl_vector*, void*);

int stls_f_reshaped (const gsl_vector*, void*, gsl_vector*);
int stls_df_reshaped (const gsl_vector*, void*, gsl_matrix*);
int stls_fdf_reshaped (const gsl_vector*, 
	      void*, gsl_vector*, gsl_matrix*);


void cholesky_of_block_of_reshaped_gamma( stls_opt_data_reshaped* );
void jacobian_reshaped( stls_opt_data_reshaped*,  gsl_matrix*);



/* Old functions */
void allocate_and_prepare_data_old( gsl_matrix* a, gsl_matrix* b, const data_struct* s, opt_and_info *opt, stls_opt_data_old *P );
void free_memory_old( stls_opt_data_old *P );


void cholgam( stls_opt_data_old* );
void jacobian( stls_opt_data_old*,  gsl_matrix*);
void xmat2xext( gsl_matrix_const_view, gsl_matrix *, int);


double stls_f_ (const gsl_vector*, void*);

int stls_f (const gsl_vector*, void*, gsl_vector*);
int stls_df (const gsl_vector*, void*, gsl_matrix*);
int stls_fdf (const gsl_vector*, 
	      void*, gsl_vector*, gsl_matrix*);


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



