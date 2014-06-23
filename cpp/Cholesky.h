/** Abstract class for Cholesky factorization of \f$\mathrm{L}_{\Gamma}^T \mathrm{L}_{\Gamma} = \Gamma(R)\f$.
 * Once created, the object allocates necessary memory for computing \f$\mathrm{L}_{\Gamma} \f$ 
 * for a given Structure object. The method Cholesky::calcGammaCholesky() computes \f$\mathrm{L}_{\Gamma}\f$
 * for a given \f$R\f$ an stores it inside the object. Operations with \f$\mathrm{L}_{\Gamma}\f$ can be performed
 * using methods Cholesky::multInvCholeskyVector and Cholesky::multInvGammaVector.
 */
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