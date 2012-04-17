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
  int ls_correction;  /* Use correction computation in nonlinear least squares */
  int gcd;         /* Is gcd being computed? */

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
#define SLRA_DEF_ls_correction 0
#define SLRA_DEF_gcd          0

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
            slraAssignDefOptValue(opt, ls_correction); \
            slraAssignDefOptValue(opt, gcd); \
          } while(0)
          
#define mymax(a, b) ((a) > (b) ? (a) : (b)) 
#define mymin(a, b) ((a) < (b) ? (a) : (b))


         

#include "slra_basic.h"
#include "slra_striped.h"

#include "slra_layered_hankel.h"
#include "slra_cholesky_btbanded.h"
#include "slra_dgamma_btbanded.h"


#include "slra_func.h"


#include "slra_common.h"
#include "slralapack.h"


int slra_gsl_optimize( slraCostFunction *F, opt_and_info *opt, gsl_vector* x_vec, gsl_matrix *v );

/* Prototypes of functions */
int slra( const gsl_vector *p_in, slraStructure* s, int r, opt_and_info* opt,
          gsl_matrix *r_ini, gsl_matrix *perm, 
          gsl_vector *p_out, gsl_matrix *rh, gsl_matrix *vh );



#endif /* _SLRA_H_ */



