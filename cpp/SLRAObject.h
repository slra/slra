#include "slra.h"

class SLRAObject {
  Structure *myS;
  VarproFunction *myF;
  static void myErrorH( const char *reason, const char *F, int ln, int gsl_err );
  static gsl_error_handler_t *old_gsl_err_h;
  static size_t myObjCnt;
public:
  SLRAObject( gsl_vector p_in, gsl_vector ml, gsl_vector nk,
              gsl_matrix perm, gsl_vector wk, gsl_vector rvec,
              bool isgcd = false );
  virtual ~SLRAObject();
    
  Structure *getS() { return myS; }
  VarproFunction *getF() { return myF; }
    
  /** Run optimization
   * @param [in]     s             Structure specification
   * @param [in]     d             Rank reduction
   * @param [in,out] opt           OptimizationOptions object
   * @param [in]     Rini          Matrix for initial approximation
   * @param [in]     Psi           \f$\Psi\f$ matrix
   *                               (identity if <tt>Psi == NULL</tt> )
   * @param [out]    p_out         Approximation \f$\widehat{p}\f$
   *                               (not computed if <tt>p_out == NULL</tt> )
   * @param [out]    R_out         Output parameter vector
   *                               (not computed if <tt>R_out == NULL</tt> )
   * @param [out]    v_out         Covariance matrix for X
   * @param [out]    Rs            Matrix of vectorized Rs at each iteration
   *                               (not computed if <tt>Rs == NULL</tt> )
   * @param [out]    info          Matrix of info (time, fmin, ...)
   *                               (not computed if <tt>info == NULL</tt> )
   */
  void optimize( OptimizationOptions* opt, gsl_matrix *Rini, gsl_matrix *Psi,
             gsl_vector *p_out, gsl_matrix *r_out, gsl_matrix *v_out,
             gsl_matrix *Rs = NULL, gsl_matrix *info = NULL );
  
};

