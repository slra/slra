/* size of the work array for mb02gd */
#define EITER 1 /* maximum number of iterations reached */

/** @memberof OptimizationOptions 
 * @name Output options
 * @{*/
#define SLRA_OPT_DISP_OFF      0  /**< Disable display */
#define SLRA_OPT_DISP_FINAL    1  /**< Show only final information  */
#define SLRA_OPT_DISP_NOTIFY   2  /**< Show final information and mesages */
#define SLRA_OPT_DISP_ITER     3  /**< Display all */
/* @}*/

/** @memberof OptimizationOptions 
 * @name Optimization methods (and submethods) from GSL library.
 * For description of these methods see relevant sections of
 * <a href="http://www.gnu.org/software/gsl/manual/">GSL manual</a>
 *@{*/
/** Nonlinear Least-Squares Fitting (gsl_multifit_fdf_solver_...).
 * If opt.method == SLRA_OPT_METHOD_LM,
 * the Levenberg-Marquardt method implemented in the GSL library is used, see
 * <a href="https://www.gnu.org/software/gsl/manual/html_node/Minimization-Algorithms-using-Derivatives.html">
 * 38.8 Minimization Algorithms using Derivatives</a> in the GSL documentation.
 *
 * The algorithm is determined by the value of opt.submethod (SLRA_OPT_SUBMETHOD_QN_xxx) */
#define SLRA_OPT_METHOD_LM 0
#define SLRA_OPT_SUBMETHOD_LM_LMDER  0 /**< gsl_multifit_fdf_solver_lmder */
#define SLRA_OPT_SUBMETHOD_LM_LMSDER 1 /**< gsl_multifit_fdf_solver_lmsder */
/** Multidimensional Minimization with Derivatives (gsl_multimin_fdf_minimizer_...).
 * If opt.method == SLRA_OPT_METHOD_QN,
 * quasi-Newton methods implemented in the GSL library are used, see
 * <a href="https://www.gnu.org/software/gsl/manual/html_node/Multimin-Algorithms-with-Derivatives.html">
 * 36.7 Algorithms with Derivatives</a> in the GSL documentation.
 *
 * The algorithm is determined by the value of opt.submethod (SLRA_OPT_SUBMETHOD_QN_xxx) */
#define SLRA_OPT_METHOD_QN 1
#define SLRA_OPT_SUBMETHOD_QN_BFGS          0 /**< gsl_multimin_fdf_minimizer_bfgs */
#define SLRA_OPT_SUBMETHOD_QN_BFGS2         1 /**< gsl_multimin_fdf_minimizer_bfgs2 */
#define SLRA_OPT_SUBMETHOD_QN_CONJUGATE_PR  2 /**< gsl_multimin_fdf_minimizer_pr */
#define SLRA_OPT_SUBMETHOD_QN_CONJUGATE_FR  2 /**< gsl_multimin_fdf_minimizer_fr */
/** Multidimensional Minimization with Derivatives (gsl_multimin_fdf_minimizer_...).
 * If opt.method == SLRA_OPT_METHOD_QN,
 * Nelder-Mead method implemented in the GSL library is used, see
 * <a href="https://www.gnu.org/software/gsl/manual/gsl-ref.html#Multimin-Algorithms-without-Derivatives.html">
 * 36.8 Algorithms with Derivatives</a> in the GSL documentation.
 *
 * The algorithm is determined by the value of opt.submethod (SLRA_OPT_SUBMETHOD_NM_xxx) */
#define SLRA_OPT_METHOD_NM 2
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX       0 /**< gsl_multimin_f_minimizer_nmsimplex */
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX2      1 /**< gsl_multimin_f_minimizer_nmsimplex2 */
#define SLRA_OPT_SUBMETHOD_NM_SIMPLEX2_RAND 2 /**< gsl_multimin_f_minimizer_nmsimplex2_rand */
/** Nonlinear Least-Squares Fitting  -
 * Levenberg-Marquardt method (own implementation using the pseudoinverse).
 *
 * This is an implementation of \cite paduart10 (page 137), except that the step
 * \f$ \lambda = \frac{1}{2} \lambda\f$ is replaced by \f$ \lambda = 0.4 \lambda\f$.
 * The method uses SVD for calculation of pseudoinverse and is able to handle
 * rank-deficient Jacobians (in the case of overparameterization).
 *
 * By default, the Jacobian is scaled (normalized), as suggested in \cite paduart10.
 * The unscaled version is available is the submethod 
 * SLRA_OPT_SUBMETHOD_LMPINV_UNSCALED is selected.
 */
#define SLRA_OPT_METHOD_LMPINV 3
/**
 * At each iteration, the Jacobian is normalizes (columns are normalized).
 * This is analogous to \ref SLRA_OPT_SUBMETHOD_LM_LMSDER.
 */
#define SLRA_OPT_SUBMETHOD_LMPINV_SCALED   0
/**
 * At each iteration, the Jacobian is not scaled.
 * This is analogous to \ref SLRA_OPT_SUBMETHOD_LM_LMDER.
 */
#define SLRA_OPT_SUBMETHOD_LMPINV_UNSCALED 1

/*@}*/
 
/** @memberof OptimizationOptions 
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
#define SLRA_DEF_epscov   1e-5
#define SLRA_DEF_reggamma 0.000
#define SLRA_DEF_ls_correction 0
#define SLRA_DEF_avoid_xi 0
/* @} */


class IterationLogger {
public:
  virtual void reportIteration( int no, const gsl_vector *x, double fmin, 
                                   const gsl_vector *grad ) = 0;
};

/** Optimization options structure.
 * This structure contains input and output parameters 
 * for 'gsl_optimize' function. */
class OptimizationOptions {
public:
  /** Constructor */
  OptimizationOptions();

  /** Main function that runs GSL optimization
   * @param [in]     F     Nonlinear least squares function
   * @param [in,out] x_vec Vector containing initial approximation and returning
   *                       the minimum point 
   * @param [out]    v     Covariance matrix for x
   */
  int gslOptimize( NLSFunction *F, gsl_vector* x_vec, gsl_matrix *v,
                   IterationLogger *itLog );

  /** Main function that runs LM optimization (for the method SLRA_OPT_METHOD_LMPINV)
   * @param [in]     F     Nonlinear least squares function
   * @param [in,out] x_vec Vector containing initial approximation and returning
   *                       the minimum point 
   */
  int lmpinvOptimize( NLSFunction *F, gsl_vector* x_vec, IterationLogger *itLog );

  /** Initialize method and submethod fields from string 
   * @param [in]     str   a string consisting of one or two characters
   *                       
   * The first character determines the value of opt.method.
   * | str[0] |  value of opt.method
   * |--------|-----------------------------
   * |   'l'  | \ref SLRA_OPT_METHOD_LM
   * |   'q'  | \ref SLRA_OPT_METHOD_QN
   * |   'n'  | \ref SLRA_OPT_METHOD_NM
   * |   'p'  | \ref SLRA_OPT_METHOD_LMPINV
   *
   * The second determines the value of opt.submethod:
   * | str[0] | str[1] | value of opt.submethod
   * |--------|--------|---------------------------------
   * |   'l'  | 'l'    | \ref SLRA_OPT_SUBMETHOD_LM_LMDER
   * |   'l'  | 's'    | \ref SLRA_OPT_SUBMETHOD_LM_LMSDER
   * |   'q'  | 'b'    | \ref SLRA_OPT_SUBMETHOD_QN_BFGS
   * |   'q'  | '2'    | \ref SLRA_OPT_SUBMETHOD_QN_BFGS2
   * |   'q'  | 'p'    | \ref SLRA_OPT_SUBMETHOD_QN_CONJUGATE_PR
   * |   'q'  | 'f'    | \ref SLRA_OPT_SUBMETHOD_QN_CONJUGATE_FR
   * |   'n'  | 'n'    | \ref SLRA_OPT_SUBMETHOD_NM_SIMPLEX
   * |   'n'  | '2'    | \ref SLRA_OPT_SUBMETHOD_NM_SIMPLEX2
   * |   'n'  | 'r'    | \ref SLRA_OPT_SUBMETHOD_NM_SIMPLEX2_RAND
   * |   'p'  | 's'    | \ref SLRA_OPT_SUBMETHOD_LMPINV_SCALED
   * |   'p'  | 'u'    | \ref SLRA_OPT_SUBMETHOD_LMPINV_UNSCALED
   * if the second letter is absent the first submethod is selected.
   */
  void str2Method( const char *str );

  /** @name General-purpose options  */  
  ///@{
  int disp;  //!< displayed information, see SLRA_OPT_DISP_xxx
  int method;    ///< method, see SLRA_OPT_METHOD_xxx 
  int submethod; ///< submethod, see SLRA_OPT_SUBMETHOD_<method>_xxx 
  ///@}
  
  /** @name Stopping criteria parameters */  
  ///@{
  size_t maxiter;///< Maximal number of iterations
  double epsabs; ///< epsabs in gsl_multifit_test_delta (see GSL documentation)
  double epsrel; ///< epsrel in gsl_multifit_test_delta (see GSL documentation)
  double epsgrad;///< epsabs in gsl_multimin_test_gradient or 'gsl_multifit_test_gradient'
  double epsx;   ///< epsabs in gsl_multimin_test_size  (used only in Nelder-Mead)
  double maxx;   ///< Maximum absolute value of the elements of the parameter vector
  ///@}
  
  /** @name Method-specific parameters */  
  ///@{
  double step;   ///< 'step_size' for fdfminimizer_set, fminimizer_set 
  double tol;    ///< 'tol' for fdfminimizer_set, fminimizer_set
  double epscov; ///< Eps for cutoff when computing covariance matrix
  ///@}
  
  /** @name Advanced parameters */  
  ///@{
  double reggamma;   ///< regularization parameter for gamma, absolute 
  int ls_correction; ///< Use correction computation in Levenberg-Marquardt 
  int avoid_xi;      ///< Avoid [X I] representation, and use own Levenberg-Marquardt
  ///@}

  /** @name Output info */  
  ///@{
  size_t iter;  ///< Total number of iterations
  double fmin;  ///< Value of the cost function 
  double time;  ///< Time spent on local optimization 
  ///@}
};






