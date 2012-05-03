class Exception {
  static const int MSG_MAX = 200;

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
  virtual void calcGammaCholesky( gsl_matrix *R ) = 0;
  /** Solve linear system with factor \f$C\f$.
   * Computes \f$ y_r \leftarrow C^{-1} y_r\f$  if trans = 0 or 
   * \f$ y_r \leftarrow C^{-T} y_r\f$  if trans = 1 */
  virtual void multInvCholeskyVector( gsl_vector * yr, int trans ) = 0;  
  /** Solves linear system with factor \f$C\f$. 
   * Computes \f$M_r^T \leftarrow C^{-1} M_r^T\f$  if trans = 0 or 
   * \f$M_r^T \leftarrow C^{-T} M_r^T\f$  if trans = 1 */
  virtual void multInvCholeskyTransMatrix( gsl_matrix * M_r, int trans );
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
  virtual void calcYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, 
                   gsl_vector *yr ) = 0;
  /** Calculate the nonlinear part of pseudojacobian.
   * Calculates \f$\displaystyle res \leftarrow \frac{\partial}{\partial X_{ij}} 
    \Gamma \left(\Phi \begin{bmatrix}X\\-I_d\end{bmatrix} \right)y_r\f$
   * */
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr ) = 0;
};


/** Abstract class for problem structure specification.
 * This class includes information about matrix structure and weights.
 */
class Structure {
public:
  virtual ~Structure() {}
  virtual int getNp() const = 0;     ///< Returns \f$n_p\f$
  virtual int getNplusD() const = 0; ///< Returns \f$m\f$
  virtual int getM() const = 0;      ///< Returns \f$n\f$
  
  /** Fills matrix from given parameter vector */
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  = 0; 
  /** Corrects a given vector p.
   * Computes \f$ p \leftarrow p - G^T(R) y_r\f$,
   * where \f$y_r\f$ is a precomputed \f$\Gamma^{-1}(R) s(R)\f$. */
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr,
                         bool scaled = true ) = 0;
  /** Creates Cholesky object for this structure and rank reduction */
  virtual Cholesky *createCholesky( int D, double reg_gamma ) const = 0;
  /** Creates DGamma object for this structure and rank reduction */
  virtual DGamma *createDGamma( int D ) const = 0;
};

/** Abstract class for s-dependent structure.
 * Assumes that the rows of the matrix are s-dependent, i.e.
 * \f${\bf cov}(S(\widetilde{p})_i, S(\widetilde{p})_j) =const\cdot W_{i,j}\f$
 * and \f$W_{i,j} = 0\f$ for \f$|i-j| < s\f$,
 * where \f$\wtilde{p}\f$ is an uncorrelated sequence with 
 * \f${\bf D} \widetilde{p}_k = \gamma_k\f$.
 *
 * In particular, this implies that \f$\Gamma(R)\f$ is a block matrix with
 * blocks \f$R^{\rm T} W_{i,j} R \f$.
 */
class SDependentStructure : public Structure {
public:
  /** Returns \f$s\f$. */
  virtual int getS() const = 0; 
  
  /** Returns \f$res \leftarrow W_{i,j} B\f$ */
  virtual void WijB( gsl_matrix *res, int i, int j, 
                     const gsl_matrix *B ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{i,j} B\f$ */
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{i,j} V\f$ */
  virtual void AtWijV( gsl_vector *res, int i, int j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const = 0;
};


/** Abstract class for stationary s-dependent structure.
 * A subclass of s-dependent structure, where
 * \f$W_{i,j} = W_{j-i}\f$.
 *
 * In other words, sequence of rows \f$S(\widetilde{p})_j\f$ is
 * stationary.
 */
class StationaryStructure : public SDependentStructure {

public:
  /** Returns \f$res \leftarrow W_{k} B\f$ */
  virtual void WkB( gsl_matrix *res, int k, const gsl_matrix *B ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{k} B\f$ */
  virtual void AtWkB( gsl_matrix *res, int k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWkB, double beta = 0 ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{k} B\f$ */
  virtual void AtWkV( gsl_vector *res, int k,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWkV, double beta = 0 ) const = 0;
                      
  virtual void WijB( gsl_matrix *res, int i, int j, 
                     const gsl_matrix *B ) const {
    WkB(res, j- i, B);
  }
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWijB, double beta = 0 ) const {
    AtWkB(res, j - i, A, B, tmpWijB, beta);
  }

  virtual void AtWijV( gsl_vector *res, int i, int j, 
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const {
    AtWkV(res, j - i, A, V, tmpWijV, beta);
  }                      
};
