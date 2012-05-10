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


/** @memberof OptimizationOptions 
 * @name Output options
 * @{*/
#define SLRA_OPT_DISP_NOTIFY   0  /**< Display information mesages */
#define SLRA_OPT_DISP_FINAL    1  /**< Show only final information ? */
#define SLRA_OPT_DISP_ITER     2  /**< Show each iteration */
#define SLRA_OPT_DISP_OFF      3  /**< Disable display */
/* @}*/

/** @memberof OptimizationOptions 
 * @name Optimization methods (and submethods) from GSL library.
 * For description of these methods see relevant sections of
 * <a href="http://www.gnu.org/software/gsl/manual/">GSL manual</a>
 *@{*/
/** Nonlinear Least-Squares Fitting (gsl_multifit_fdf_solver_...) -
 * Levenberg-Marquardt method */
#define SLRA_OPT_METHOD_LM 0
#define SLRA_OPT_SUBMETHOD_LM_LMDER  0 /**< ..._lmder */
#define SLRA_OPT_SUBMETHOD_LM_LMSDER 1 /**< ..._lmsder */
/** Multidimensional Minimization with Derivatives 
 * (gsl_multifit_fdf_solver_...) -
 * quasi-Newton and conjugate gradients methods */
#define SLRA_OPT_METHOD_QN 1
#define SLRA_OPT_SUBMETHOD_QN_BFGS          0 /**< ..._bfgs */
#define SLRA_OPT_SUBMETHOD_QN_BFGS2         1 /**< ..._bfgs2 */
#define SLRA_OPT_SUBMETHOD_QN_CONJUGATE_PR  2 /**< ..._pr */
#define SLRA_OPT_SUBMETHOD_QN_CONJUGATE_FR  2 /**< ..._fr */
/** Multidimensional Minimization without Derivatives 
 * (gsl_multifit_fdf_solver_...) -
 * Nead-Melder method */
#define SLRA_OPT_METHOD_NM 2
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX       0 /**< ..._nmsimplex */
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX2      1 /**< ..._nmsimplex2 */
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX2_RAND 2 /**< ..._nmsimplex2_rand */
/*@}*/
 
/** @memberof OptimizationOptions 
 * @defgroup t1
 * @name Default values for parameters
 * @{ */
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
/* @} */

/** Optimization options structure.
 * This structure contains input and output parameters 
 * for 'gsl_optimize' function. */
class OptimizationOptions {
public:
  /** Constructor */
  OptimizationOptions() : disp(SLRA_DEF_disp), method(SLRA_DEF_method),
      submethod(SLRA_DEF_submethod),  maxiter(SLRA_DEF_maxiter),
      epsabs(SLRA_DEF_epsabs), epsrel(SLRA_DEF_epsrel), 
      epsgrad(SLRA_DEF_epsgrad), epsx(SLRA_DEF_epsx),
      step(SLRA_DEF_step), tol(SLRA_DEF_tol), reggamma(SLRA_DEF_reggamma),
      ls_correction(SLRA_DEF_ls_correction), gcd(SLRA_DEF_gcd) {
  }
  
  /** Initialize disp field from string */
  void str2Disp( const char *str );
  /** Initialize method and submethod fields from string */
  void str2Method( const char *str );

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
};

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
#define PSI_STR "Psi"
          
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


/** Main function that runs SLRA optimization
 * @ingroup MainFunctions 
 * @param [in]  p_in  Input pararmeter vector \f$p\f$  
 * @param [in]  s     Structure specification  
 * @param [in]  r     Structure specification  
 * @param [in]  opt   Optimization options
 * @param [in]  Rini  Initial approximation  
 * @param [in]  Phi   Phi matrix
 *                    (identity if Phi == NULL )   
 * @param [out] p_out Approximation \f$\widehat{p}\f$  
 *                    (not computed if p_out == NULL )
 * @param [out] R_out Output parameter vector 
 */
int slra( const gsl_vector *p_in, Structure* s, int r, 
          OptimizationOptions* opt, gsl_matrix *Rini, gsl_matrix *Phi, 
          gsl_matrix *Psi, gsl_vector *p_out, gsl_matrix *rh, gsl_matrix *vh );
/** Function that runs SLRA optimization */
int gsl_optimize( CostFunction *F, OptimizationOptions *opt, 
                       gsl_vector* x_vec, gsl_matrix *v );
/** @defgroup MainFunctions
 * Global functions. */


#endif /* _SLRA_H_ */



