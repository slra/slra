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

#elif defined(BUILD_MEX_OCTAVE) || defined(BUILD_MEX_MATLAB)

#include "mex.h"
#define PRINTF mexPrintf
#define WARNING mexWarnMsgTxt

#else

#include <stdio.h>
#define PRINTF printf
#define WARNING printf

#endif

#define DEBUGINT(x) PRINTF("%s = %d\n", #x, x)
#define DEBUGDOUBLE(x) PRINTF("%s = %f\n", #x, x)

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

extern char meth_codes[];
extern char *submeth_codes[];

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
  double inv_w;            /* Square root of inverse of the weight */
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



class slraGammaComputations {
public:  
  virtual  ~slraGammaComputations() {}
  virtual void computeCholeskyOfGamma( gsl_matrix *R ) = 0;

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans ) = 0;  
  virtual void multiplyInvGammaVector( gsl_vector * yr ) = 0;                
  virtual void multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans ) { /* Default implementation */
    for (int i = 0; i < yr_matr->size1; i++) {
      gsl_vector_view row = gsl_matrix_row(yr_matr, i);
      multiplyInvCholeskyVector(&row.vector, trans);
    }
  }
};

class slraDerivativeComputations {
public:  
  virtual ~slraDerivativeComputations() {}
  virtual void computeYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, gsl_vector *yr ) = 0;

  virtual void computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr ) = 0;

};


class slraStructure {
public:
  virtual ~slraStructure() {}
  virtual int getNp() const = 0;
  virtual int getNplusD() const = 0;
  virtual int getM() const = 0;
  virtual void fillMatrixFromP( gsl_matrix* c, gsl_vector* p )  = 0; 
  
  virtual slraGammaComputations *createGammaComputations( int r, double reg_gamma ) = 0;
  virtual slraDerivativeComputations *createDerivativeComputations( int r ) = 0;
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) = 0;
};



class slraWkInterface {
public:
  virtual const gsl_matrix *getWk( int k )  const = 0; 
  virtual int getS() const  = 0;
  virtual int getNplusD() const = 0;
};




class slraFlexStructure : virtual public slraStructure, virtual public slraWkInterface {
  int myK;                      /* = rowdim(block in T/H blocks) */ 
  int myQ;	                /* number of blocks in C = [C1 ... Cq] */
  
  int myNp;
  
  int myNplusD;
  int myNpScale, myNpOffset;
  int myMaxLag;

  void computeStats();
  void computeWkParams(); 

  gsl_matrix **myA;
  

  slraFlexBlock *mySA;	/* q-element array describing C1,...,Cq; */  
public:
  slraFlexStructure( const slraFlexStructure &s ); /* Copy constructor */
  slraFlexStructure( const data_struct *s, int np = -1 ); /* "Copy" constructor */
  slraFlexStructure( const double *s_matr, int q, int k, int s_matr_cols, int np_or_m = -1, bool set_m = false, 
                     const double *w_k = NULL );
  virtual ~slraFlexStructure();

  virtual int getNp() const { return myNp; }
  virtual int getNplusD() const { return myNplusD; }
  virtual int getM() const { return (myNp - myNpOffset) / myNpScale; }
  virtual slraGammaComputations *createGammaComputations( int r, double reg_gamma );
  virtual slraDerivativeComputations *createDerivativeComputations( int r );

  void setNp( int np );
  void setM( int m );
  
  int getQ() const { return myQ; }
  int getK() const { return myK; }

  int getMaxLag() const { return myMaxLag; }

  int getNpOffset() const { return myNpOffset; }
  int getNpScale() const { return myNpScale; }

//  const slraFlexBlock & getFlexBlock( int i ) const { return mySA[i]; }
  
  int getFlexBlockLag( int l ) const { return mySA[l].blocks_in_row; }
  int getFlexBlockNCol( int l ) const { return mySA[l].blocks_in_row * mySA[l].nb; }
  int getFlexBlockNb( int l ) const { return mySA[l].nb; }
  bool isFlexBlockExact( int l ) const { return mySA[l].exact; }
  bool isFlexBlockToeplitz( int l ) const { return mySA[l].toeplitz; }
  double getInvBlockWeight( int l ) const { return mySA[l].inv_w; }
  
  int getFlexBlockT( int l ) const { return getFlexBlockLag(l) + (getM() / getK()) - 1; }
  int getFlexBlockNp( int l ) const { return getFlexBlockT(l) * getK() * getFlexBlockNb(l); }
  
  virtual void fillMatrixFromP( gsl_matrix* c, gsl_vector* p ); 
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );


  virtual int getS() const { return myMaxLag; }
  virtual const gsl_matrix *getWk( int k ) const { 
    return myA[k]; 
  }
};


class slraFlexGammaComputations : virtual public slraGammaComputations {
protected:
  int my_use_slicot;
  double my_reg_gamma;
  
  const slraWkInterface *myW;
  
  int myN, myD;
  
  int myMg;
  
  int s_minus_1;
  int d_times_s;
  int d_times_Mg;
  int d_times_s_minus_1;
  
  int myCholeskyWorkSize;
 
  double *myGammaVec;
  gsl_matrix *myGamma;
  gsl_matrix *myWkTmp;
  double *myPackedCholesky;
  double *myCholeskyWork;
  
public:
  slraFlexGammaComputations( const slraWkInterface *s, int r, int Mg,
     int use_slicot, double reg_gamma  );
  virtual ~slraFlexGammaComputations();

  int getD() const { return myD; }

  virtual void computeCholeskyOfGamma( gsl_matrix *R );

  virtual void multiplyInvPartCholeskyArray( double * yr, int trans, int size, int chol_size );
  virtual void multiplyInvPartGammaArray( double * yr, int size, int chol_size );


  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
    if (yr->stride != 1) {
      throw new slraException("Cannot multiply vectors with stride != 1\n");
    }
    multiplyInvPartCholeskyArray(yr->data, trans, yr->size, d_times_Mg);
  }
  virtual void multiplyInvGammaVector( gsl_vector * yr ) {
    if (yr->stride != 1) {
      throw new slraException("Cannot multiply vectors with stride != 1\n");
    }
    multiplyInvPartGammaArray(yr->data, yr->size, d_times_Mg);
  }

  virtual void multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans ) {
    if (yr_matr->size2 != yr_matr->tda) {
      slraGammaComputations::multiplyInvCholeskyTransMatrix(yr_matr, trans);
    } else {
      multiplyInvPartCholeskyArray(yr_matr->data, trans, yr_matr->size1 * yr_matr->size2, d_times_Mg);
    }
  }

};

class slraFlexDerivativeComputations : virtual public 
                                       slraDerivativeComputations {
  const slraWkInterface *myW;
  int  myD, myK;
  
  gsl_vector *myTempWkColRow;
  gsl_matrix *myDGamma;
  
  gsl_matrix *myWk_R;
  gsl_matrix *myWkT_R;
  gsl_matrix *myN_k;

public:
  slraFlexDerivativeComputations( const slraWkInterface *s, int r, int K );
  virtual ~slraFlexDerivativeComputations();

  int getD() const { return myD; }

  
  virtual void computeYrtDgammaYr( gsl_matrix *mgrad_r, gsl_matrix *R, gsl_vector *yr );
  virtual void computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );

};


class slraFlexStructureExt : public slraStructure {
  slraFlexStructure mySimpleStruct;

  int myN;
  int *myOldMl;
//  double *myWk;
  int myMaxMl;
  
  
  /* Helper variables */
  int myM;
  int myNp;

public:
  slraFlexStructureExt( int q, int N, double *oldNk, double *oldMl, double *Wk );
  virtual ~slraFlexStructureExt();
  
  virtual int getM() const { return myM; }
  virtual int getNp() const { 
    return mySimpleStruct.getNpOffset() * myN + mySimpleStruct.getNpScale() * myM; 
  }
  virtual int getNplusD() const { return mySimpleStruct.getNplusD(); }


  virtual slraGammaComputations *createGammaComputations( int r, double reg_gamma );
  virtual slraDerivativeComputations *createDerivativeComputations( int r );


  int getMl( int k ) const { return myOldMl[k]; }
  int getMaxMl() const { return myMaxMl; }
  int getBlocksN() const { return myN; }

  
  const slraWkInterface *getWkInterface() const { 
    return &mySimpleStruct; 
  }

  virtual void fillMatrixFromP( gsl_matrix* c, gsl_vector* p ) ;
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );
};



class slraFlexGammaComputationsExt : public slraGammaComputations {
  slraFlexGammaComputations myBase;
  slraFlexStructureExt *myStruct;
public:  
  slraFlexGammaComputationsExt( slraFlexStructureExt *s, int r, int use_slicot, double reg_gamma  ) :
      myStruct(s), myBase(s->getWkInterface(), r, s->getMaxMl(), use_slicot, reg_gamma) {
  }
  virtual ~slraFlexGammaComputationsExt() {}
  
  virtual void computeCholeskyOfGamma( gsl_matrix *R ) {
    myBase.computeCholeskyOfGamma(R);
  }
  
  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multiplyInvGammaVector( gsl_vector * yr );                
  

};



class slraFlexDerivativeComputationsExt : virtual public 
                                       slraDerivativeComputations {
  slraFlexDerivativeComputations myBase;
  slraFlexStructureExt *myStruct;
  gsl_matrix *myTmpGrad;
public:  
  slraFlexDerivativeComputationsExt( slraFlexStructureExt *s, int r  ) ;
  virtual ~slraFlexDerivativeComputationsExt();

  virtual void computeYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, gsl_vector *yr );

  virtual void computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );

};




class slraCostFunction {
  slraStructure *myStruct;
  int myRank;
  slraGammaComputations *myGam;
  slraDerivativeComputations *myDeriv;
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

  slraCostFunction( slraStructure *s, int r, gsl_vector *p, opt_and_info *opt, gsl_matrix *perm  );
  virtual ~slraCostFunction();
  
  int getD() { return myStruct->getNplusD() - myRank; }
  int getNplusD() { return myStruct->getNplusD(); }
  int getN() { return myRank; }
  int getM() { return myStruct->getM(); }

  const gsl_matrix * getPerm() { return myPerm; }
  const gsl_matrix * getSMatr() { return myMatr; }
  
  
  void computeR( gsl_matrix_const_view x_mat, gsl_matrix *R ); 
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
  void computeSr( gsl_matrix *R, gsl_vector *Sr );


  void computeRGammaSr( const gsl_vector *x, gsl_matrix *R, gsl_vector *Sr ) {
    computeR(x, myTmpR);
    myGam->computeCholeskyOfGamma(myTmpR);
    computeSr(myTmpR, Sr);
  } 
  
  void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *R, gsl_matrix *jac );
  void computeGradFromYr( gsl_vector* yr, gsl_matrix *R, gsl_vector *grad );


  void computeFuncAndPseudoJacobianLs( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac );
  void computeFuncAndGrad( const gsl_vector* x, double * f, gsl_vector *grad );
  
  void  computeCorrection( gsl_vector* p, const gsl_vector* x );
  
  static int slra_f_ls( const gsl_vector* x, void* params, gsl_vector *res ) {
    ((slraCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int slra_df_ls( const gsl_vector* x,  void* params, gsl_matrix *jac ) {
    ((slraCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int slra_fdf_ls( const gsl_vector* x,  void* params, gsl_vector *res, gsl_matrix *jac ) {
    ((slraCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, jac);
    return GSL_SUCCESS;
  }


  static double  slra_f( const gsl_vector* x, void* params ) {
    double f;
    ((slraCostFunction *)params)->computeFuncAndGrad(x, &f, NULL);
    return f;
  }
  static void  slra_df( const gsl_vector* x, void* params, gsl_vector *grad ) {
    ((slraCostFunction *)params)->computeFuncAndGrad(x, NULL, grad);
  }
  static void  slra_fdf( const gsl_vector* x, void* params, double *f, gsl_vector *grad ) {
    ((slraCostFunction *)params)->computeFuncAndGrad(x, f, grad);
  }
};

int slra_gsl_optimize( slraCostFunction *F, opt_and_info *opt, gsl_vector* x_vec, gsl_matrix *v );

/* Prototypes of functions */
int slra(gsl_vector* p, slraStructure * s, int rank, gsl_matrix* x,
         gsl_matrix* v, opt_and_info* opt, int x_given, int compute_ph,
         gsl_matrix* perm);
	

void tmv_prod_new(gsl_matrix*, int, 
	      gsl_vector*, int, gsl_vector*);

int tls(gsl_matrix*, gsl_matrix*, gsl_matrix*);
 
#endif

#ifdef __cplusplus
extern "C" {
#endif

void print_mat(const gsl_matrix*);
void print_mat_tr(const gsl_matrix*);
void print_arr(double*, int);

void gsl_matrix_vectorize(double*, gsl_matrix*);
void gsl_matrix_vec_inv(gsl_matrix*, double*);


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



