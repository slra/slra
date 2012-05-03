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


#include <gsl/gsl_blas.h>

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

/** Optimization options structure.
 * This structure contains input and output parameters 
 * for 'gsl_optimize' function. */
typedef struct {
  /** @name General-purpose options  */  
  ///@{
  int disp;  //!< displayed information, see SLRA_OPT_DISP_xxx
  int method;    ///< method, see SLRA_OPT_METHOD_xxx 
  int submethod; ///< submethod, see SLRA_OPT_SUBMETHOD_<method>_xxx 
  ///@}
  
  /** @name Stopping criteria parameters */  
  ///@{
  int maxiter;           ///< Maximal number of iterations 
  double epsabs, epsrel; ///< Eps for 'gsl_multifit_test_delta' criterion
  double epsgrad;        ///< Eps for 'gsl_multi..._test_gradient' criteria
  double epsx;           ///< Eps for Nead-Melder method
  ///@}
  
  /** @name Method-specific parameters */  
  ///@{
  double step;  ///< 'step_size' for fdfminimizer_set, fminimizer_set 
  double tol;   ///< 'tol' for fdfminimizer_set, fminimizer_set
  ///@}
  
  /** @name Advanced parameters */  
  ///@{
  double reggamma;   ///< regularization parameter for gamma, absolute 
  int ls_correction; ///< Use correction computation in Levenberg-Marquardt 
  int gcd;           ///< Testing option for agcd problem 
  ///@}

  /** @name Output info */  
  ///@{
  int iter;     ///< Total number of iterations 
  double fmin;  ///< Value of the cost function 
  double time;  ///< Time spent on local optimization 
  ///@}
} OptimizationOptions;

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
#define SLRA_DEF_ls_correction 0
#define SLRA_DEF_gcd          0

#define AssignDefOptValue(opt, field) \
  do { opt.field = SLRA_DEF_##field; } while(0)

#define AssignDefOptValues(opt) do {  \
            AssignDefOptValue(opt, disp); \
            AssignDefOptValue(opt, method); \
            AssignDefOptValue(opt, submethod);   \
            AssignDefOptValue(opt, maxiter); \
            AssignDefOptValue(opt, epsabs); \
            AssignDefOptValue(opt, epsrel); \
            AssignDefOptValue(opt, epsgrad); \
            AssignDefOptValue(opt, epsx); \
            AssignDefOptValue(opt, step); \
            AssignDefOptValue(opt, tol); \
            AssignDefOptValue(opt, reggamma); \
            AssignDefOptValue(opt, ls_correction); \
            AssignDefOptValue(opt, gcd); \
          } while(0)
          
/* field names for s */
#define ML_STR "m"
#define NK_STR "n"
#define PERM_STR "phi"
#define WK_STR "w"

/* field names for opt */
#define RINI_STR "Rini"
#define DISP_STR "disp"
#define METHOD_STR "method"

/* names for output */
#define PH_STR "ph"
#define INFO_STR "info"


/* field names for info */
#define RH_STR "Rh"
#define VH_STR "Vh"
#define FMIN_STR "fmin"
#define ITER_STR "iter"
#define TIME_STR "time"
          
#define mymax(a, b) ((a) > (b) ? (a) : (b)) 
#define mymin(a, b) ((a) < (b) ? (a) : (b))

#include "slra_basic.h"
#include "slra_striped.h"

#include "slra_layered_hankel.h"
#include "slra_layered_hankel_weighted.h"
#include "slra_cholesky_btbanded.h"
#include "slra_dgamma_btbanded.h"

#include "slra_computation.h"

#include "slra_common.h"
#include "slralapack.h"


/** \addtogroup Func
 * Global functions. */
///{@

/** Main function that runs SLRA optimization */
int slra( const gsl_vector *p_in, Structure* s, int r, 
          OptimizationOptions* opt, gsl_matrix *r_ini, gsl_matrix *perm, 
          gsl_vector *p_out, gsl_matrix *rh, gsl_matrix *vh );
/** Function that runs SLRA optimization */
int gsl_optimize( CostFunction *F, OptimizationOptions *opt, 
                       gsl_vector* x_vec, gsl_matrix *v );

///@}

#endif /* _SLRA_H_ */



