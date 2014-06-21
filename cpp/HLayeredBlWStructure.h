/** Layered Hankel structure with blockwise weights.
 * The layered Hankel structure is a layered structure with \f$q\f$ Hankel blocks:
 * \f[ \mathscr{H}_{{\bf m}, n} := 
 * \begin{bmatrix}
 * \mathscr{H}_{m_1,n} (p^{(1)}) \\  \vdots \\ \mathscr{H}_{m_q,n} (p^{(q)}) 
 * \end{bmatrix}, 
 * \f]
 * see eqn. \f$\mathscr{H}_{{\bf m}, n}\f$ in \cite slra-efficient. Here 
 * \f$p = \mathrm{col}(p^{(1)},\ldots,p^{(q)})\f$.
 *
 * The blockwise weights are the weights 
 * \f[\mathrm{col}(\omega_1,\ldots,\omega_q),\f]
 * defined in eqn. \f$(\mbox{block-wise }w)\f$ \cite slra-efficient.
 */
class HLayeredBlWStructure : public StationaryStructure {
  typedef struct {
    size_t blocks_in_row;       /* Number of blocks in a row of Ci */
    double inv_w;            /* Square root of inverse of the weight */
  } Layer;

  size_t myQ;	                /* number of layers */
  size_t myN;
  size_t myM;
  size_t myMaxLag;
  gsl_matrix **myA;
  Layer *mySA;	/* q-element array describing C1,...,Cq; */  

  void computeStats();
  void computeWkParams(); 
protected:
  size_t nvGetNp() const { return (myN - 1) * myQ + myM; }  
public:
  /** Constructs the HLayeredBlWStructure object.
   * @param[in] m_vec \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}^{\top}\f$
   * @param[in] q     number of blocks \f$q\f$
   * @param[in] n     number of columns \f$n\f$
   * @param[in] w_vec vector of weights 
   * \f${\bf w} =\begin{bmatrix} \omega_1 & \cdots & \omega_q \end{bmatrix}^{\top}\f$.
   * If `w_vec == NULL` then \f${\bf w}\f$ is set to be
   * \f$\begin{bmatrix}1&\cdots&1\end{bmatrix}^{\top}\f$. */
  HLayeredBlWStructure( const double *m_vec, size_t q, size_t n, 
                        const double *w_vec = NULL );
  virtual ~HLayeredBlWStructure();
  /** @name Implementing Structure interface */
  /**@{*/
  virtual size_t getM() const { return myM; }
  virtual size_t getN() const { return myN; }
  virtual size_t getNp() const { return nvGetNp(); }
  virtual Cholesky *createCholesky( size_t d ) const;
  virtual DGamma *createDGamma( size_t d ) const;
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ); 
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
                                   const gsl_vector *y, 
                                   double alpha = -1, double beta = 1,
                                   bool skipFixedBlocks = true ); 
  virtual void multByWInv( gsl_vector* p, long deg = 2 );
  /**@}*/
 
  /** @name Implementing StationaryStructure interface */
  /**@{*/
  virtual size_t getS() const { return myMaxLag; }
  virtual void WkB( gsl_matrix *X, long k, const gsl_matrix *B ) const;
  virtual void AtWkB( gsl_matrix *X, long k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpVkB, double beta = 0 ) const;
  virtual void AtWkV( gsl_vector *u, long k,
                      const gsl_matrix *A, const gsl_vector *v, 
                      gsl_vector *tmpVkV, double beta = 0 ) const;
  /**@}*/
 
  /** @name HLayeredBlWStructure-specific methods */
  /**@{*/
  /** Returns \f$q\f$. */
  size_t getQ() const { return myQ; }
  /** Returns \f$\max \{m_l\}_{l=1}^{q}\f$ */
  size_t getMaxLag() const { return myMaxLag; }
  /** Returns \f$m_{l}\f$.
   * Internally, the indexing of blocks starts from 0. 
   * @param[in] l_1 \f$0\f$-based index, such that \f$l = l_1+1\f$ and \f$0\le l_1 < q\f$. */
  size_t getLayerLag( size_t l_1 ) const { return mySA[l_1].blocks_in_row; }
  /** Checks whether \f$w_{l}=\infty\f$ (block is fixed).
   * @copydetails HLayeredBlWStructure::getLayerLag */
  bool isLayerExact( size_t l_1 ) const { return (mySA[l_1].inv_w == 0.0); }
  /** Returns \f$w_{l}^{-1}\f$ 
   * @copydetails HLayeredBlWStructure::getLayerLag */
  double getLayerInvWeight( size_t l_1 ) const { return mySA[l_1].inv_w; }
  /** Returns number of parameters in each block: \f$n_p^{(l)} = m_l+n-1\f$ 
   * @copydetails HLayeredBlWStructure::getLayerLag */
  size_t getLayerNp( size_t l_1 ) const { return getLayerLag(l_1) + getN() - 1; }
  /** Returns a pointer to \f$\mathrm{V}_{k}\f$ matrix.  @todo remove.
   *  @param[in] k index \f$k: 0\le k \le \mu\f$. */
  const gsl_matrix *getWk( size_t k ) const { return myA[k]; } 
  /**@}*/
};


