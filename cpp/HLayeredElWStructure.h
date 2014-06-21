/** Layered Hankel structure with elementwise weights.
 * The layered Hankel structure \f$\mathscr{H}_{{\bf m}, n}\f$
 * is defined in description of HLayeredBlWStructure.
 *
 * The elementwise weights are the weights 
 * \f[\mathrm{col}(w_1,\ldots,w_{n_p}),\f]
 * corresponding to the case 3  of Theorem 1 in \cite slra-efficient.
 *
 * The computations are based on the expression of \f$\mathrm{V}_{\#ij}\f$ 
 * as in \f$(\mathrm{V}_{\#ij}(\mathscr{H}_{{\bf m}, n}))\f$ 
 */
class HLayeredElWStructure : public MuDependentStructure {
  HLayeredBlWStructure myBase;
  gsl_vector *myInvWeights;
  gsl_vector *myInvSqrtWeights;
  void mulInvWij( gsl_matrix * res, long i_1 ) const;
public:  
  /** Constructs WLayeredHStructure object.
   * @param m_vec \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}^{\top}\f$
   * @param w_vec vector of weights 
   * \f${\bf w} =\begin{bmatrix} w_1 & \cdots & w_{n_p} \end{bmatrix}^{\top}\f$.
   */
  HLayeredElWStructure( const double *m_vec, size_t q, size_t n, 
                        const double *w_vec );
  virtual ~HLayeredElWStructure();

  /** @name Implementing Structure interface */
  /**@{*/
  virtual size_t getM() const { return myBase.getM(); }
  virtual size_t getN() const { return myBase.getN(); }
  virtual size_t getNp() const { return myBase.getNp(); }
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p );
  virtual Cholesky *createCholesky( size_t d ) const;
  virtual DGamma *createDGamma( size_t d ) const;
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
                                   const gsl_vector *y, 
                                   double alpha = -1, double beta = 1,
                                   bool skipFixedBlocks = true ); 
  virtual void multByWInv( gsl_vector* p, long deg = 2 );
  /**@}*/
  
  /** @name Implementing SDependentStructure interface */
  /**@{*/
  virtual size_t getMu() const { return myBase.getMu(); }
  virtual void VijB( gsl_matrix *X, long i_1, long j_1, 
                     const gsl_matrix *B ) const;
  virtual void AtVijB( gsl_matrix *X, long i_1, long j_1, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpVijB, double beta = 0 ) const;
  virtual void AtVijV( gsl_vector *u, long i_1, long j_1,
                      const gsl_matrix *A, const gsl_vector *v, 
                      gsl_vector *tmpVijV, double beta = 0 ) const;
  /**@}*/
                      
  /** @name LayeredHStructure-specific methods */
  /**@{*/
  /** @brief @copybrief HLayeredBlWStructure::getQ() */
  size_t getQ() const { return myBase.getQ(); }
  /** @brief @copybrief HLayeredBlWStructure::getMaxLag() */
  size_t getMaxLag() const { return myBase.getMaxLag(); }
  /** @brief @copybrief HLayeredBlWStructure::getLayerLag() 
   * @copydetails HLayeredBlWStructure::getLayerLag */
  size_t getLayerLag( size_t l_1 ) const { return myBase.getLayerLag(l_1); }
  /** @brief @copybrief HLayeredBlWStructure::getLayerNp()
   * @copydetails HLayeredBlWStructure::getLayerLag */
  size_t getLayerNp( size_t l_1 ) const { return myBase.getLayerNp(l_1); }
  /** Returns \f$w_i\f$ 
   * @param[in] i_1 \f$0\f$-based index, such that \f$i = i_1+1\f$ and 
   * \f$0\le i_1 < n_p\f$. */
  double getInvWeight( size_t i_1 ) const { 
    return gsl_vector_get(myInvWeights, i_1); 
  }
  /**@}*/
};






