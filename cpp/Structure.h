/** Abstract class for specification of matrix structure and weights.
 * This class contains information about:
 *  - matrix structure \f$\mathscr{S}: \mathbb{R}^{n_p} \to \mathbb{R}^{m \times n}\f$,  
 *    (eq. \f$(\mathscr{S})\f$) in  \cite slra-efficient );
 *  - weight matrix \f$\mathrm{W}: \mathbb{R}^{n_p \times n_p}\f$, 
 *     (eq. \f$(\|\cdot\|^2_{\mathrm{W}})\f$ in \cite slra-efficient ).
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
  
  /** Updates \f$p\f$ as \f$ p \leftarrow \beta p + \alpha \left((I_n \otimes R) 
   ** \mathbf{S}_{\mathscr{S}}\right) y\f$.
   * The matrix \f$(I_n \otimes R) \mathbf{S}_{\mathscr{S}}\f$ is the unweighted 
   * part of \f$G(R)\f$ defined in eq. (\f$G(R)\f$) in \cite slra-efficient.
   * \param[out,in]  p   parameter vector \f$p\in\mathbb{R}^{n_p}\f$
   * \param[in]  Rt      matrix \f$R^{\top} \in \mathbb{R}^{d \times m} \f$,
   *                     where \f$d\f$d  is the rank reduction.  
   * \param[in]  y       vector \f$y \in \mathbb{R}^{nd} \f$ 
   * \param[in]  alpha   constant \f$\alpha\f$ 
   * \param[in]  beta   constant \f$\beta\f$ 
   * \param[in]  skipFixedBlocks   @todo: explain 
   */
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
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

