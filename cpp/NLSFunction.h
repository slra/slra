/** Abstract class for nonlinear least squares functions
 * Represents a function of the form
 * \f$f(x) = \|g(x)\|_2^2\f$, where \f$g: \mathbb{R}^{n} \to \mathbb{R}^{n_s}\f$,
 * i.e., a sum of squares of univariate functions.
 */
class NLSFunction {
public:
  virtual ~NLSFunction() {}
  /** Returns \f$n\f$ --- the number of variables */
  virtual size_t getNvar() = 0;
  /** Returns the number of essential variables */
  virtual size_t getNEssVar() { getNvar(); }
  /** Returns \f$n_s\f$ --- the number of squares (dimension of \f$g\f$) */
  virtual size_t getNsq() = 0;
  /** Computes function \f$f\f$ and gradient */
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, 
                                   gsl_vector *grad ) = 0;
  /** Computes the vector \f$g\f$ and the Jacobian (or pseudo-jacobian) */
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, 
                                  gsl_matrix *jac ) = 0;

  static double _f( const gsl_vector* x, void* params ) {
    double f;
    ((NLSFunction *)params)->computeFuncAndGrad(x, &f, NULL);
    return f;
  }
  static void _df( const gsl_vector* x, void* params, gsl_vector *grad ) {
    ((NLSFunction *)params)->computeFuncAndGrad(x, NULL, grad);
  }
  static void _fdf( const gsl_vector* x, void* params, double* f, 
                     gsl_vector* grad ) {
    ((NLSFunction *)params)->computeFuncAndGrad(x, f, grad);
  }

  static int _f_ls( const gsl_vector* x, void* params, gsl_vector* res ) {
    ((NLSFunction *)params)->computeFuncAndJac(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int _df_ls( const gsl_vector* x,  void* params, gsl_matrix* jac ) {
    ((NLSFunction *)params)->computeFuncAndJac(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int _fdf_ls( const gsl_vector* x, void* params, gsl_vector*res, 
                       gsl_matrix *jac ) {
    ((NLSFunction *)params)->computeFuncAndJac(x, res, jac);
    return GSL_SUCCESS;
  }
};



class NLSVarpro : public NLSFunction {
public: 
  /** Returns \f$n\f$ --- the number of variables */
  virtual size_t getNvar() = 0;
  /** Returns \f$n_s\f$ --- the number of squares (dimension of \f$g\f$) */
  virtual size_t getNsq() = 0;
    
  /** Computes function \f$f\f$ and gradient */
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, 
                                   gsl_vector *grad ) = 0;
  /** Computes the vector \f$g\f$ and the Jacobian (or pseudo-jacobian) */
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, 
                                  gsl_matrix *jac ) = 0;

  virtual size_t getD() = 0;
  virtual size_t getM() = 0;
								  
  virtual void RTheta2x( gsl_matrix *RTheta, gsl_vector *x ) = 0; 
  virtual void x2RTheta( gsl_matrix *RTheta, const gsl_vector *x ) = 0; 
  
  /* To remove: */
  virtual void computePhat( gsl_vector* p, const gsl_vector* x ) = 0;
  virtual void computeDefaultx( gsl_vector *x ) = 0; 

  virtual gsl_matrix x2xmat( const gsl_vector *x ) = 0;
};