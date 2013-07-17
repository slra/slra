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
  /** Computes \f$ p \leftarrow \beta p + \alpha G^T(R) y\f$. 
   * The fixed blocks may be skipped. */
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *R, 
                                   const gsl_vector *y, double alpha = -1,
                                   double beta = 0,
                                   bool skipFixedBlocks = true ) = 0; 
  /** Computes \f$ p \leftarrow diag(w)^{-deg/2} p \f$. */
  virtual void multByWInv( gsl_vector* p, long deg = 2 ) = 0;
  /** Creates Cholesky object for this structure and rank reduction */
  virtual Cholesky *createCholesky( size_t D ) const = 0;
  /** Creates DGamma object for this structure and rank reduction */
  virtual DGamma *createDGamma( size_t D ) const = 0;
};

