class Exception {
  static const size_t MSG_MAX = 200;

  char myMsg[MSG_MAX];
public:
  Exception( const char *msg, ... );
  const char *getMessage()  { return myMsg; }
};


/** Abstract class for Cholesky factorization of \f$\Gamma(R)\f$.
 * Computation Cholesky factorization \f$C^T C = \Gamma(R)\f$ */
class Cholesky {
public:  
  virtual  ~Cholesky() {}
  
  /** Computes Cholesky factorization */
  virtual void calcGammaCholesky( const gsl_matrix *R, double reg_gamma = 0 ) = 0;
  /** Solve linear system with factor \f$C\f$.
   * Computes \f$ y_r \leftarrow C^{-1} y_r\f$  if trans = 0 or 
   * \f$ y_r \leftarrow C^{-T} y_r\f$  if trans = 1 */
  virtual void multInvCholeskyVector( gsl_vector * yr, long trans ) = 0;  
  /** Solves linear system with factor \f$C\f$. 
   * Computes \f$M_r^T \leftarrow C^{-1} M_r^T\f$  if trans = 0 or 
   * \f$M_r^T \leftarrow C^{-T} M_r^T\f$  if trans = 1 */
//  virtual void multInvCholeskyTransMatrix( gsl_matrix * M_r, long trans );
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


/** Abstract class for problem structure specification.
 * This class includes information about matrix structure and weights.
 */
class Structure {
public:
  virtual ~Structure() {}
  virtual size_t getNp() const = 0;     ///< Returns \f$n_p\f$
  virtual size_t getM() const = 0;      ///< Returns \f$m\f$
  virtual size_t getN() const = 0;      ///< Returns \f$n\f$
  
  /** Fills matrix from given parameter vector */
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  = 0; 
  /** Corrects a given vector p.
   * Computes \f$ p \leftarrow p - G^T(R) y_r\f$,
   * where \f$y_r\f$ is a precomputed \f$\Gamma^{-1}(R) s(R)\f$. */
  virtual void correctP( gsl_vector* p, const gsl_matrix *R, const gsl_vector *yr,
                         long wdeg = 2 ) = 0;
  /** Creates Cholesky object for this structure and rank reduction */
  virtual Cholesky *createCholesky( size_t D ) const = 0;
  /** Creates DGamma object for this structure and rank reduction */
  virtual DGamma *createDGamma( size_t D ) const = 0;
};




