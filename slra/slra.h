/*********************************
 * Header file fo SLRA package
 *********************************/

/* slra.h: SLRA header file */
#ifndef _SLRA_H_
#define _SLRA_H_

#if defined(BUILD_R_PACKAGE)

#include <R.h>
#define PRINTF Rprintf
#define WARNING Rprintf

#elif defined(BUILD_MEX_OCTAVE) ||  defined(BUILD_MEX_MATLAB)

#include "mex.h"
#define PRINTF mexPrintf
#define WARNING mexWarnMsgTxt

#else

#include <stdio.h>
#define PRINTF printf
#define WARNING printf

#endif

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h> /* Levenberge-Marquardt */
#include <gsl/gsl_multimin.h>      /* BFGS Newton-type     */

/* size of the work array for mb02gd */
#define EITER 1 /* maximum number of iterations reached */

#ifdef __cplusplus
extern "C" {
#endif


#define SLRA_OPT_DISP_NOTIFY   0
#define SLRA_OPT_DISP_FINAL    1
#define SLRA_OPT_DISP_ITER     2
#define SLRA_OPT_DISP_OFF      3



#define SLRA_OPT_METHOD_LM   0
#define SLRA_OPT_METHOD_QN   1
#define SLRA_OPT_METHOD_NM   2

#define SLRA_OPT_SUBMETHOD_LM_LMDER         0
#define SLRA_OPT_SUBMETHOD_LM_LMSDER        1


#define SLRA_OPT_SUBMETHOD_QN_BFGS          0
#define SLRA_OPT_SUBMETHOD_QN_BFGS2         1
#define SLRA_OPT_SUBMETHOD_QN_CONJUGATE_PR  2
#define SLRA_OPT_SUBMETHOD_QN_CONJUGATE_FR  3

#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX       0
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX2      1
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX2_RAND 2

/* optimization options and output information structure */
typedef struct {
  int disp; /* displayed information: 
	       1 - notify, 2 - final, 3 - iter, 4 - off */
  
  /* method */
  int method;
  int submethod;
  
  /* stopping criterion */  
  int maxiter;
  double epsabs, epsrel;    /* Eps for LM */
  double epsgrad;           /* Eps for LM and QN */
  double epsx;              /* Eps for NM */ 
  
  /* optimization parameters */
  double step;
  double tol;
  
  double reggamma; /* To be worked out */
  int use_slicot;

  /* output information */
  int iter;
  double fmin;
  double time;
  double chol_time;
} opt_and_info;

#define SLRA_DEF_disp       SLRA_OPT_DISP_NOTIFY 
#define SLRA_DEF_method     SLRA_OPT_METHOD_LM
#define SLRA_DEF_submethod  0
#define SLRA_DEF_maxiter  100 
#define SLRA_DEF_epsabs   0
#define SLRA_DEF_epsrel   1e-5
#define SLRA_DEF_epsgrad  1e-5
#define SLRA_DEF_epsx     1e-5
#define SLRA_DEF_step     0.001
#define SLRA_DEF_tol      1e-6
#define SLRA_DEF_reggamma 0.001
#define SLRA_DEF_use_slicot 1

#define slraAssignDefOptValue(opt, field) \
  do { opt.field = SLRA_DEF_##field; } while(0)

#define slraAssignDefOptValues(opt) do {  \
            slraAssignDefOptValue(opt, disp); \
            slraAssignDefOptValue(opt, method); \
            slraAssignDefOptValue(opt, submethod);   \
            slraAssignDefOptValue(opt, maxiter); \
            slraAssignDefOptValue(opt, epsabs); \
            slraAssignDefOptValue(opt, epsrel); \
            slraAssignDefOptValue(opt, epsgrad); \
            slraAssignDefOptValue(opt, epsx); \
            slraAssignDefOptValue(opt, step); \
            slraAssignDefOptValue(opt, tol); \
            slraAssignDefOptValue(opt, reggamma); \
            slraAssignDefOptValue(opt, use_slicot); \
          } while(0)
          
          

/* structure in the data matrix C = [ A B ] */ 
#define MAXQ 10	/* maximum number of blocks in C */

typedef struct {
  int blocks_in_row;       /* Number of blocks in a row of Ci */
  int nb;                  /* Number of columns in each small block */
  int exact;               /* 1 - exact block, 0 - not exact */  
  int toeplitz;            /* 1 - Toeplitz marix, 0 - Hankel matrix */  
} slraFlexBlock;

typedef struct {
  int k;	          /* = rowdim(block in T/H blocks) */ 
  int q;	          /* number of blocks in C = [C1 ... Cq] */
  slraFlexBlock a[MAXQ];  /* q-element array describing C1,...,Cq; */  
} data_struct;

/* additional info about matrix structure */
typedef struct {
  int total_cols;
  int np_scale, np_offset; /* Coefficients in the formula  
			      np = np_scale * m + np_offset */
} flex_struct_add_info;

/* three dimensional array of covariances w */
typedef struct {
  int s;	/* length of the array (w.s = s+1 from the paper) */
  gsl_matrix **a;
} w_data;

#define mymax(a, b) ((a) > (b) ? (a) : (b)) 
#define mymin(a, b) ((a) < (b) ? (a) : (b))

#define COMMON_PARAMS  \
    int m, n, d; \
    gsl_matrix* c; \
    double reggamma; \
    int k;  \
    int  n_plus_d, 		/* = col_dim(C) */ \
    n_times_d,			/* = number of elements in x */ \
    k_times_d,			/* = row_dim(gamma) */ \
    k_times_d_times_s,		/* = col_dim(gamma) */ \
    k_times_d_times_s_minus_1,  /* = col_dim(gamma) - 1 */ \
    m_times_d, 			/* = row_dim(rb) */ \
    m_div_k,  \
    s_minus_1;  \
    int one; /* One for blas routines */ \
    int use_slicot; \
    int chol_count; \
    clock_t chol_time; 

#define PREPARE_COMMON_PARAMS(C, Nn, S, OPT, PP, isblock) \
  do {\
  PP->m = C->size1; \
  PP->n = Nn; \
  PP->d = C->size2 - Nn;		\
  /* set other parameters */ \
  PP->c = gsl_matrix_alloc(C->size1, C->size2); \
  gsl_matrix_memcpy(PP->c, C); \
  PP->reggamma = OPT->reggamma; \
  /* find Wk  */ \
  s2w(S, &PP->w, PP->n+PP->d, isblock);\
  PP->k = S->k;  \
  PP->n_plus_d = PP->n + PP->d;   \
  PP->n_times_d = PP->n * PP->d;   \
  PP->k_times_d = S->k * PP->d; \
  PP->k_times_d_times_s = PP->k_times_d * PP->w.s;\
  PP->k_times_d_times_s_minus_1 = PP->k_times_d_times_s - 1;\
  PP->m_times_d = PP->m * PP->d;\
  PP->m_div_k = (int) PP->m / S->k;\
  PP->s_minus_1 = PP->w.s - 1;\
  PP->one = 1;\
  PP->use_slicot = OPT->use_slicot; \
  PP->chol_count = 0; \
  PP->chol_time = 0; \
  } while (0)

typedef struct {
  COMMON_PARAMS;

  w_data w; \

  /* Preallocated arrays */  
  gsl_matrix *x_ext; 
  gsl_vector *yr;
  
  /* Preallocated arrays for cholgam */
  double *rb;   /* Result of Cholesky factorization */
  gsl_matrix *tmp; /* Temp matrix for cholgam (x_ext' * w_k) 
		      P->k_times_d x SIZE_W  */
  gsl_matrix *gamma;
  double *gamma_vec;
  int ldwork;       /* Size of Dwork for MB02GD  */
  double *dwork;    /* Dwork for MB02GD  */

  /* Preallocated arrays for jacobian */
  gsl_matrix *dgamma, *st;
  double *jres1, * jres2;
} slra_opt_data_old;


#ifdef __cplusplus
}


class slraException {
  static const int MSG_MAX = 200;

  char myMsg[MSG_MAX];
public:
  slraException( const char *msg, ... );
  const char *getMessage()  { return myMsg; }
};




class slraFlexStructure {
  int myK;                      /* = rowdim(block in T/H blocks) */ 
  int myQ;	                /* number of blocks in C = [C1 ... Cq] */
  
  int myNp;
  
  int myNplusD;
  int myNpScale, myNpOffset;
  int myMaxLag;
  
  void computeStats();
  

  slraFlexBlock mySA[MAXQ];	/* q-element array describing C1,...,Cq; */  
public:
  slraFlexStructure( const data_struct *s, int np = -1 ); /* "Copy" constructor */
  slraFlexStructure( const double *s_matr, int q, int k, int s_matr_cols, int np = -1 );
  virtual ~slraFlexStructure() {}

  void setNp( int np, bool check_solvability = true );
  int getQ() const { return myQ; }
  int getK() const { return myK; }
  int getNp() const { return myNp; }
  int getNplusD() const { return myNplusD; }
  int getM() const { return (myNp - myNpOffset) / myNpScale; }
  int getMaxLag() const { return myMaxLag; }
  
//  const slraFlexBlock & getFlexBlock( int i ) const { return mySA[i]; }
  
  int getFlexBlockLag( int l ) const { return mySA[l].blocks_in_row; }
  int getFlexBlockNCol( int l ) const { return mySA[l].blocks_in_row * mySA[l].nb; }
  int getFlexBlockNb( int l ) const { return mySA[l].nb; }
  bool isFlexBlockExact( int l ) const { return mySA[l].exact; }
  bool isFlexBlockToeplitz( int l ) const { return mySA[l].toeplitz; }
  
  int getFlexBlockT( int l ) const { return getFlexBlockLag(l) + (getM() / getK()) - 1; }
  int getFlexBlockNp( int l ) const { return getFlexBlockT(l) * getK() * getFlexBlockNb(l); }
  
  void fillMatrixFromP( gsl_matrix* c, gsl_vector* p ) const ; 
  

};


class slraGammaComputations {
public:  
  virtual  ~slraGammaComputations() {}
  virtual const double *getPackedCholesky() = 0;
  virtual void computeCholeskyOfGamma( gsl_matrix *R ) = 0;
  virtual void multiplyInvCholeskyArray( double * yr, int trans, int rep = 1 ) = 0;  /*  yr - I/O */
 
  virtual void multiplyInvGammaArray( double * yr ) = 0;  /* yr - I/O */

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
    multiplyInvCholeskyArray(yr->data, trans, 1);
  }
  virtual void multiplyInvGammaVector( gsl_vector * yr ){
    multiplyInvGammaArray(yr->data);
  }
  
  virtual void correctVector( gsl_vector* p, slraFlexStructure *s, 
                                               gsl_matrix *R, gsl_vector *f ) = 0;
};


class slraFlexComputationsParams {
  int myS;	/* length of the array (w.s = s+1 from the paper) */
  gsl_matrix **myA;

public:
  slraFlexComputationsParams( const slraFlexStructure *s ); 
  virtual ~slraFlexComputationsParams(); 
  
  const int getS() { return myS; }
  const gsl_matrix *getWk( int k ) { return myA[k]; }
};



class slraDerivativeComputations {
public:  
  virtual ~slraDerivativeComputations() {}
  virtual void computeYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, gsl_vector *yr ) = 0;

  virtual void computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr ) = 0;

};



class slraFlexGammaComputations : virtual public slraGammaComputations {
protected:
  int my_use_slicot;
  double my_reg_gamma;
  
  slraFlexComputationsParams *myW;
  int myM, myN, myD, myK;
  
  int m_div_k;
  int s_minus_1;
  int d_times_s;
  int d_times_m_div_k;
  int d_times_s_minus_1;
  
  int myCholeskyWorkSize;
 
  double *myGammaVec;
  gsl_matrix *myGamma;
  gsl_matrix *myWkTmp;
  double *myPackedCholesky;
  double *myCholeskyWork;
  
public:
  slraFlexGammaComputations( slraFlexStructure *s, int r, 
     int use_slicot, double reg_gamma, slraFlexComputationsParams *w  );
  virtual ~slraFlexGammaComputations();

  virtual const double *getPackedCholesky() { return myPackedCholesky; }
  virtual void computeCholeskyOfGamma( gsl_matrix *R );
  virtual void multiplyInvCholeskyArray( double * yr, int trans, int rep = 1 );
  virtual void multiplyInvGammaArray( double * yr );
  virtual void correctVector( gsl_vector* p, slraFlexStructure *s, 
                                               gsl_matrix *R, gsl_vector *f );
};


class slraFlexDerivativeComputations : virtual public 
                                       slraDerivativeComputations {
  slraFlexComputationsParams *myW;
  int myM, myN, myD, myK;
  
  int d_times_m_div_k;
  int m_div_k;
  
  gsl_vector *myTempWkColRow;
  gsl_matrix *myDGamma;
  
  gsl_matrix *myWk_R;
  gsl_matrix *myWkT_R;
  gsl_matrix *mySubPhiTmp;
  gsl_matrix *myN_k;

public:
  slraFlexDerivativeComputations( slraFlexStructure *s, int r, 
      slraFlexComputationsParams *w  );
  virtual ~slraFlexDerivativeComputations();
  
  virtual void computeYrtDgammaYr( gsl_matrix *mgrad_r, gsl_matrix *R, gsl_vector *yr );
  virtual void computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );

};



class slraFlexCostFunction;


/* data needed for cost function and Jacobian evaluation */
typedef struct {
  COMMON_PARAMS;

  /* Preallocated arrays for cholgam (new) */
  int  d_times_s;		/* = col_dim(new gamma) */
  int  d_times_s_minus_1;  /* = col_dim(new gamma) - 1 */
  int  d_times_m_div_k;  /* */

  double *rb2;   /* Result of Cholesky factorization (test) */

  /* Preallocated data for block computations */
  gsl_matrix *bx_ext;  

  double *brg_rb;   /* Result of Cholesky factorization */
  gsl_matrix *brg_tmp; /* Temp matrix for cholgam (x_ext' * w_k) 
			  P->k_times_d x SIZE_W  */
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
  
  gsl_matrix *brg_c;
  
  double *brg_j1b_vec;
  gsl_matrix *brg_j1b;
  gsl_matrix *brg_j1b_2;

  gsl_vector *brg_j1_cvec;
  gsl_vector *brg_j2_pvec;

  
  /* Helper functions for gradient */
  gsl_matrix *brg_grad_N_k;
  gsl_matrix *brg_grad_Vx_k;
  gsl_matrix *brg_grad_tmp1;
  gsl_matrix *brg_grad_tmp2;
  gsl_matrix *brg_grad_tmp3;

  gsl_matrix *perm;  
  
  gsl_matrix *brg_perm_tmp;  
  
  slraFlexComputationsParams *myW;
  slraGammaComputations *myGamma;
  slraDerivativeComputations *myDerivative;
  slraFlexStructure *myStruct;
  
  
  slraFlexCostFunction *myCostFun;
} slra_opt_data_reshaped;


class slraFlexCostFunction {
//  slraFlexStructure * const myStruct;
  slraFlexStructure myStructure;
  int myRank;
  slraFlexComputationsParams myW;
  slraFlexGammaComputations myGamma;
  slraFlexDerivativeComputations myDerivative;
  gsl_matrix *myMatr;
  gsl_matrix *myPerm;
  
  gsl_matrix *myMatrMulPerm;
  
  /* Helper computation variables */
  gsl_matrix *myTmpGradR;
  gsl_matrix *myTmpGradR2;
  gsl_matrix *myTmpR;  
  gsl_vector *myTmpYr;  

  /* Jacobian computation */
  double *myTmpJacobianArray;  
  gsl_vector *myTmpJacobianCol;  

  gsl_matrix *myTmpGrad;  

public:

  slraFlexCostFunction( slraFlexStructure s, int r, 
      gsl_vector *p, opt_and_info *opt, gsl_matrix *perm  );
  virtual ~slraFlexCostFunction();
  
  int getD() { return myStructure.getNplusD() - myRank; }
  int getNplusD() { return myStructure.getNplusD(); }
  int getN() { return myRank; }
  int getM() { return myStructure.getM(); }

  const gsl_matrix * getPerm() { return myPerm; }
  const gsl_matrix * getSMatr() { return myMatr; }
  
  
  void computeR( gsl_matrix_const_view x_mat, gsl_matrix *R ); 
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
  void computeSr( gsl_matrix *R, gsl_vector *Sr );


  void computeRGammaSr( const gsl_vector *x, gsl_matrix *R, gsl_vector *Sr ) {
    computeR(x, myTmpR);
    myGamma.computeCholeskyOfGamma(myTmpR);
    computeSr(myTmpR, Sr);
  } 
  
  void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *R, gsl_matrix *jac );
  void computeGradFromYr( gsl_vector* yr, gsl_matrix *R, gsl_vector *grad );


  void computeFuncAndPseudoJacobianLs( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac );
  void computeFuncAndGrad( const gsl_vector* x, double * f, gsl_vector *grad );
  
  void  computeCorrection( gsl_vector* p, const gsl_vector* x );
  
  static int slra_f_ls( const gsl_vector* x, void* params, gsl_vector *res ) {
    ((slraFlexCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int slra_df_ls( const gsl_vector* x,  void* params, gsl_matrix *jac ) {
    ((slraFlexCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int slra_fdf_ls( const gsl_vector* x,  void* params, gsl_vector *res, gsl_matrix *jac ) {
    ((slraFlexCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, jac);
    return GSL_SUCCESS;
  }


  static double  slra_f( const gsl_vector* x, void* params ) {
    double f;
    ((slraFlexCostFunction *)params)->computeFuncAndGrad(x, &f, NULL);
    return f;
  }
  static void  slra_df( const gsl_vector* x, void* params, gsl_vector *grad ) {
    ((slraFlexCostFunction *)params)->computeFuncAndGrad(x, NULL, grad);
  }
  static void  slra_fdf( const gsl_vector* x, void* params, double *f, gsl_vector *grad ) {
    ((slraFlexCostFunction *)params)->computeFuncAndGrad(x, f, grad);
  }
};



extern "C" {

/* TODO: replace with something that uses printf. 
   #define print_vec(v) gsl_vector_fprintf(stdout,v,"%16.14f") */

/*void xmat2_block_of_xext( gsl_matrix_const_view, gsl_matrix *,
			  gsl_matrix *, gsl_matrix *);

void allocate_and_prepare_data_reshaped( gsl_matrix* c, int n,
         opt_and_info *opt, slra_opt_data_reshaped *P, gsl_matrix *perm);
void free_memory_reshaped( slra_opt_data_reshaped *P );

double slra_f_reshaped_ (const gsl_vector*, void*);

int slra_f_reshaped (const gsl_vector*, void*, gsl_vector*);
int slra_df_reshaped (const gsl_vector*, void*, gsl_matrix*);
int slra_fdf_reshaped (const gsl_vector*, 
	      void*, gsl_vector*, gsl_matrix*);

void cholesky_of_block_of_reshaped_gamma( slra_opt_data_reshaped* );
void jacobian_reshaped( slra_opt_data_reshaped*,  gsl_matrix*);

int slra_allocate_params( void *pparams, gsl_vector* p, data_struct* s, gsl_matrix* x,
         gsl_matrix* v, opt_and_info* opt, int x_given, int compute_ph,
         gsl_matrix *perm, int perm_given  );
int slra_gsl_optimize( slra_opt_data_reshaped *P, opt_and_info *opt, gsl_vector* x_vec, gsl_matrix *v );  

void grad_reshaped( slra_opt_data_reshaped* P, gsl_vector* grad );*/

#endif


/* Prototypes of functions */
int slra(gsl_vector* p, data_struct* s, int rank, gsl_matrix* x,
         gsl_matrix* v, opt_and_info* opt, int x_given, int compute_ph,
         gsl_matrix* perm);
	
int check_and_adjust_parameters(data_struct *s, flex_struct_add_info *psi);
int slra_fill_matrix_from_p(gsl_matrix* c,  data_struct *s, gsl_vector* p);
int slra_correction_reshaped(gsl_vector* p, data_struct *s, void* params,
			     const gsl_vector* x);

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
 


/* Old functions */
void allocate_and_prepare_data_old( gsl_matrix* c, int n,
				    const data_struct* s, 
         opt_and_info *opt, slra_opt_data_old *P );
void free_memory_old( slra_opt_data_old *P );

void cholgam( slra_opt_data_old* );
void jacobian( slra_opt_data_old*,  gsl_matrix*);
void xmat2xext( gsl_matrix_const_view, gsl_matrix *, int);

double slra_f_ (const gsl_vector*, void*);

int slra_f (const gsl_vector*, void*, gsl_vector*);
int slra_df (const gsl_vector*, void*, gsl_matrix*);
int slra_fdf (const gsl_vector*, 
	      void*, gsl_vector*, gsl_matrix*);

void slra_df_reshaped_ (const gsl_vector* x, void* params,
			gsl_vector* grad);
void slra_fdf_reshaped_ (const gsl_vector* x, void* params, double *f,
			 gsl_vector* grad);



       



/* SLICOT and LAPACK functions */

void mb02gd_(char*, char*, int*, int*, int*, const int*, int*,
double*, int*, double*, int*, double*, const int*, int*);
void mb02md_(char*, int*, int*, int*, const int*, double*, int*,
double*, double*, int*, double*, int*, double*, int*, const int*, int*);
void dtbtrs_(char*, const char*, char*, int*, int*, const int*, const double*,
int*, double*, int*, int*);
void dpbtrs_(char*, int*, int*, const int*, const double*, int*, double*,
int*, int*);
void dpbtrf_(char*, int *, int *, double *, int *, int *);


void m_to_gsl_matrix(gsl_matrix* a_gsl, double* a_m);
void gsl_to_m_matrix(double* a_m, gsl_matrix* a_gsl); 

/* Convert double matrix to structure */
int slraMatrix2Struct( data_struct *s, double *s_matr, 
                       int q, int s_matr_cols );
void slraString2Method( const char *str_buf, opt_and_info *popt );
int slraString2Disp( const char *str_value );


#ifdef __cplusplus
}



#endif

#endif /* _SLRA_H_ */



