/** Abstract class for specification of matrix structure and weights.
 * This class contains information about:
 *  - matrix structure \f$\mathscr{S}: \mathbb{R}^{n_p} \to \mathbb{R}^{m \times n}\f$,  
 *    (equation \f$(\mathscr{S})\f$) in  in \cite slra-efficient );
 *  - weight matrix \f$\mathrm{W}: \mathbb{R}^{n_p \times n_p}\f$, 
 *     (equation \f$(\|\cdot\|^2_{\mathrm{W}})\f$ in \cite slra-efficient ).
 */
class Structure {
public:
  virtual ~Structure() {}
  virtual size_t getNp() const = 0;     ///< Returns \f$n_p\f$
  virtual size_t getM() const = 0;      ///< Returns \f$m\f$
  virtual size_t getN() const = 0;      ///< Returns \f$n\f$
  
  /** Creates a transpose of structured matrix from the parameter vector. 
   * \param[out] c  matrix \f$c \leftarrow \mathscr{S}^{\top}(p)\f$
   * \param[in]  p   parameter vector \f$p\in\mathbb{R}^{n_p}\f$
   */
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) = 0; 
  
  /** Modifes a vector as \f$ p \leftarrow \beta p + \alpha G^{\top}(R) y\f$. 
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

