/** Layered Hankel structure with elementwise weights.
 * @copydetails LayeredHStructure
 */
class WLayeredHStructure : public SDependentStructure {
  LayeredHStructure myBase;
  gsl_vector *myInvWeights;
  gsl_vector *myInvSqrtWeights;
  void mulInvWij( gsl_matrix * res, int i ) const;
public:  
  /** Constructs WLayeredHStructure object.
   * @param m_l \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}\f$
   * @param w vector of weights 
   * \f${\bf w} =\begin{bmatrix} w_1 & \cdots & w_{n_p} \end{bmatrix}\f$.
   * If w_l == NULL then \f${\bf w}\f$ is set to be
   * \f$\begin{bmatrix}1&\cdots&1\end{bmatrix}\f$. 
   */
  WLayeredHStructure( const double *m_l, size_t q, int n, 
                     const gsl_vector *w = NULL );
  virtual ~WLayeredHStructure();

  /** @name Implementing Structure interface */
  /**@{*/
  virtual int getM() const { return myBase.getM(); }
  virtual int getN() const { return myBase.getN(); }
  virtual int getNp() const { return myBase.getNp(); }
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) {
    myBase.fillMatrixFromP(c, p); 
  }
  virtual Cholesky *createCholesky( int D, double reg_gamma ) const;
  virtual DGamma *createDGamma( int D ) const;
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr,
                         int wdeg = 2 );
  /**@}*/
  
  /** @name Implementing SDependentStructure interface */
  /**@{*/
  virtual int getS() const { return myBase.getS(); }
  virtual void WijB( gsl_matrix *res, int i, int j, 
                     const gsl_matrix *B ) const;
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const;
  virtual void AtWijV( gsl_vector *res, int i, int j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const;
  /**@}*/
                      
  /** @name Structure-specific methods */
  /**@{*/
  /** @brief @copybrief LayeredHStructure::getQ()*/
  int getQ() const { return myBase.getQ(); }
  /** @brief @copybrief LayeredHStructure::getMaxLag()*/
  int getMaxLag() const { return myBase.getMaxLag(); }
  /** @brief @copybrief LayeredHStructure::getLayerLag()*/
  int getLayerLag( int l ) const { return myBase.getLayerLag(l); }
  /** @brief @copybrief LayeredHStructure::getLayerNp()*/
  int getLayerNp( int l ) const { return myBase.getLayerNp(l); }
  /** Returns \f$w_k\f$ */
  double getInvWeight( int i ) const { 
    return gsl_vector_get(myInvWeights, i); 
  }
  /**@}*/
};


/** Mosaic Hankel structure with elementwise weights.
 * @copydetails MosaicHStructure
 */
class WMosaicHStructure : public StripedStructure {
protected:
  static Structure **allocStripe( gsl_vector *m_l, gsl_vector *n_k,  
                gsl_vector *w );
public:
  /** Constructs WMosaicHStructure object.
   * @param m_l \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}\f$
   * @param n_k \f${\bf n} = \begin{bmatrix}n_1 & \cdots & n_N\end{bmatrix}\f$
   * @param w vector of weights 
   * \f${\bf w} =\begin{bmatrix} w_1 & \cdots & w_{n_p} \end{bmatrix}\f$.
   * If w_l == NULL then \f${\bf w}\f$ is set to be
   * \f$\begin{bmatrix}1&\cdots&1\end{bmatrix}\f$. 
   */
  WMosaicHStructure( gsl_vector *m_l, gsl_vector *n_k, gsl_vector *w );
  virtual ~WMosaicHStructure() {}
};





