/* size of the work array for mb02gd */
#define EITER 1 /* maximum number of iterations reached */

/** @memberof OptimizationOptions 
 * @name Output options
 * @{*/
#define SLRA_OPT_DISP_OFF      0  /**< Disable display */
#define SLRA_OPT_DISP_FINAL    1  /**< Show only final information ? */
#define SLRA_OPT_DISP_NOTIFY   2  /**< Show final information and mesages */
#define SLRA_OPT_DISP_ITER     3  /**< Display all */
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
 * @defgroup DefaultOptimizationOptions
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
#define SLRA_DEF_maxx     0
#define SLRA_DEF_step     0.001
#define SLRA_DEF_tol      1e-6
#define SLRA_DEF_reggamma 0.000
#define SLRA_DEF_ls_correction 0
/* @} */

/** Optimization options structure.
 * This structure contains input and output parameters 
 * for 'gsl_optimize' function. */
class OptimizationOptions {
public:
  /** Constructor */
  OptimizationOptions();

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
  size_t maxiter;           ///< Maximal number of iterations 
  double epsabs, epsrel; ///< Eps for 'gsl_multifit_test_delta' criterion
  double epsgrad;        ///< Eps for 'gsl_multi..._test_gradient' criteria
  double epsx;           ///< Eps for Nead-Melder method
  double maxx;           ///< Maximum X value
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
  ///@}

  /** @name Output info */  
  ///@{
  size_t iter;     ///< Total number of iterations 
  double fmin;  ///< Value of the cost function 
  double time;  ///< Time spent on local optimization 
  ///@}
};






