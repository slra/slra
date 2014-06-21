class Exception {
  static const size_t MSG_MAX = 200;

  char myMsg[MSG_MAX];
public:
  Exception( const char *msg, ... );
  const char *getMessage()  { return myMsg; }
};


/** Abstract class for Cholesky factorization of \f$\Gamma(R)\f$.
 * Cholesky factorization \f$\mathrm{L}_{\Gamma}^T \mathrm{L}_{\Gamma} = \Gamma(R)\f$ */
class Cholesky {
public:  
  virtual  ~Cholesky() {}
  
  /** Computes Cholesky factorization */
  virtual void calcGammaCholesky( const gsl_matrix *R, double reg_gamma = 0 ) = 0;
  /** Solve linear system with factor \f$C\f$.
   * Computes \f$ y_r \leftarrow C^{-1} y_r\f$  if trans = 0 or 
   * \f$ y_r \leftarrow C^{-T} y_r\f$  if trans = 1 */
  virtual void multInvCholeskyVector( gsl_vector * yr, long trans ) = 0;  
  /** Solves linear system with \f$\Gamma(R)\f$ . 
   * Computes \f$y_r \leftarrow \Gamma^{-1} y_r\f$ 
   * using Cholesky factorization */
  virtual void multInvGammaVector( gsl_vector * yr ) = 0;                
};

/** Abstract class for differentiating \f$\Gamma(R)\f$.
 * Computations with \f$\Gamma(R)\f$ */
class DGamma {
public:  
  virtual ~DGamma() {}
  /** Calculate the nonlinear part of the gradient.
   * Updates the gradient \f$grad \leftarrow grad + A\f$, where
   * \f$A\f$ is defined by \f$trace(A,H) = y_r d \Gamma(R, H) y_r\f$
   * */
  virtual void calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                   const gsl_vector *yr ) = 0;
  /** Calculate the nonlinear part of pseudojacobian.
   * Calculates \f$\displaystyle res \leftarrow \frac{\partial}{\partial R_{ij}} 
    \Gamma \left(R \right)y_r\f$
   * */
  virtual void calcDijGammaYr( gsl_vector *res, const gsl_matrix *R, 
                   size_t i, size_t j, const gsl_vector *Yr ) = 0;
};






