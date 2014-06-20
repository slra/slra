/** Layered Hankel structure with elementwise weights.
 * @copydetails HLayeredBlWStructure
 */
class HLayeredElWStructure : public SDependentStructure {
  HLayeredBlWStructure myBase;
  gsl_vector *myInvWeights;
  gsl_vector *myInvSqrtWeights;
  void mulInvWij( gsl_matrix * res, long i ) const;
public:  
  /** Constructs WLayeredHStructure object.
   * @param m_l \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}\f$
   * @param w vector of weights 
   * \f${\bf w} =\begin{bmatrix} w_1 & \cdots & w_{n_p} \end{bmatrix}\f$.
   */
  HLayeredElWStructure( const double *m_l, size_t q, size_t n, 
                      const double *w );
  virtual ~HLayeredElWStructure();

  /** @name Implementing Structure interface */
  /**@{*/
  virtual size_t getM() const { return myBase.getM(); }
  virtual size_t getN() const { return myBase.getN(); }
  virtual size_t getNp() const { return myBase.getNp(); }
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p );
  virtual Cholesky *createCholesky( size_t D ) const;
  virtual DGamma *createDGamma( size_t D ) const;
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
                                   const gsl_vector *y, 
                                   double alpha = -1, double beta = 1,
                                   bool skipFixedBlocks = true ); 
  virtual void multByWInv( gsl_vector* p, long deg = 2 );
  /**@}*/
  
  /** @name Implementing SDependentStructure interface */
  /**@{*/
  virtual size_t getS() const { return myBase.getS(); }
  virtual void WijB( gsl_matrix *res, long i, long j, 
                     const gsl_matrix *B ) const;
  virtual void AtWijB( gsl_matrix *res, long i, long j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const;
  virtual void AtWijV( gsl_vector *res, long i, long j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const;
  /**@}*/
                      
  /** @name Structure-specific methods */
  /**@{*/
  /** @brief @copybrief LayeredHStructure::getQ()*/
  size_t getQ() const { return myBase.getQ(); }
  /** @brief @copybrief LayeredHStructure::getMaxLag()*/
  size_t getMaxLag() const { return myBase.getMaxLag(); }
  /** @brief @copybrief LayeredHStructure::getLayerLag()*/
  size_t getLayerLag( size_t l ) const { return myBase.getLayerLag(l); }
  /** @brief @copybrief LayeredHStructure::getLayerNp()*/
  size_t getLayerNp( size_t l ) const { return myBase.getLayerNp(l); }
  /** Returns \f$w_k\f$ */
  double getInvWeight( size_t i ) const { 
    return gsl_vector_get(myInvWeights, i); 
  }
  /**@}*/
};






