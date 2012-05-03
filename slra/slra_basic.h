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
  
  //! Computes Cholesky factorization
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
  /** Calculate varying part of the gradient.
   * */
  virtual void calcYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, 
                   gsl_vector *yr ) = 0;
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr ) = 0;
};


/** Abstract class for structure specification */
class Structure {
public:
  virtual ~Structure() {}
  virtual int getNp() const = 0;     ///< Returns \f$n_p\f$
  virtual int getNplusD() const = 0; ///< Returns \f$m\f$
  virtual int getM() const = 0;      ///< Returns \f$n\f$
  
  /** Fills matrix from given parameter vector */
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  = 0; 
  /** Computes correction \f$ p \leftarrow p - G^T(R) y_r\f$ */
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr,
                         bool scaled = true ) = 0;
  
  /** Creates Cholesky object for this structure and rank reduction */
  virtual Cholesky *createCholesky( int D, double reg_gamma ) const = 0;
  /** Creates DGamma object for this structure and rank reduction */
  virtual DGamma *createDGamma( int D ) const = 0;
};

class SDependentStructure : public Structure {
public:
  virtual int getS() const = 0;
  virtual void WijB( gsl_matrix *res, int i, int j, 
                     const gsl_matrix *B ) const = 0;
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const = 0;
  virtual void AtWijV( gsl_vector *res, int i, int j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const = 0;
};

class StationaryStructure : public SDependentStructure {

public:
  virtual void WkB( gsl_matrix *res, int k, const gsl_matrix *B ) const = 0;
  virtual void AtWkB( gsl_matrix *res, int k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWkB, double beta = 0 ) const = 0;
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
