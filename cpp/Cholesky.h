/** Abstract class for Cholesky factorization of \f$\mathrm{L}_{\Gamma}^T \mathrm{L}_{\Gamma} = \Gamma(R)\f$.
 * Once created, the object allocates necessary memory for computing \f$\mathrm{L}_{\Gamma} \f$ 
 * for a given Structure object. The method Cholesky::calcGammaCholesky() computes \f$\mathrm{L}_{\Gamma}\f$
 * for a given \f$R\f$ an stores it inside the object. Operations with \f$\mathrm{L}_{\Gamma}\f$ can be performed
 * using methods Cholesky::multInvCholeskyVector and Cholesky::multInvGammaVector.
 */
class Cholesky {
public:  
  virtual  ~Cholesky() {}

  /** Computes and stores the Cholesky factor \f$\mathrm{L}_{\Gamma}\f$ of \f$\Gamma(R)\f$. 
   * If \f$\Gamma(R)\f$ is singular and \f$\gamma > 0\f$, the function tries to compute the
   * Cholesky factorization of \f$\Gamma(R)+\gamma I_{nd}\f$. If the computations fail, then
   * the function aborts all computations (throws an Exception).
   * @param[in] Rt        the matrix \f$R^{\top} \in \mathbb{R}^{m \times d}\f$.
   * @param[in] reg  a regularization parameter \f$\gamma\f$. 
   */
  virtual void calcGammaCholesky( const gsl_matrix *Rt, double reg = 0 ) = 0;

  /** Solve linear system with factor \f$\mathrm{L}_{\Gamma}\f$ (or its submatrix).
   * Computes 
   * * \f$ y_r \leftarrow (\mathrm{L}_{\Gamma}^{-1})_{1:P,1:P} y_r\f$  if `trans == 0` or 
   * * \f$ y_r \leftarrow (\mathrm{L}_{\Gamma}^{-T})_{1:P,1:P} y_r\f$  if `trans == 1`. 
   * @param[in,out] yr vector \f$y_r \in \mathbb{R}^{P}\f$, where \f${P} \le nd\f$
   * @param[in]   trans `0` or `1`
   */
  virtual void multInvCholeskyVector( gsl_vector * y_r, long trans ) = 0;  

  /** Solves linear system with \f$\Gamma(R)\f$ or (its submatrix). 
   * Computes \f$y_r \leftarrow (\Gamma_{1:P,1:P})^{-1} y_r\f$ 
   * using the stored Cholesky factor. 
   * @param[in,out] yr vector \f$y_r \in \mathbb{R}^{P}\f$, where \f${P} \le nd\f$
   */
  virtual void multInvGammaVector( gsl_vector * y_r ) = 0;                
};