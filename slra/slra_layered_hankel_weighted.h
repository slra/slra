/** Layered Hankel structure with elementwise weights.
 * A structure of form
 * \f$ \mathcal{H}_{{\bf m}, n} := 
 * \begin{bmatrix} \mathcal{H}_{n,m_1} & \cdots & \mathcal{H}_{n,m_q} 
 * \end{bmatrix} \f$  with \f$q\f$ blocks
 * and elementwise weights 
 * \f$\begin{bmatrix} w_1 & \cdots & w_{n_p} \end{bmatrix}\f$ 
 */
class WLayeredHStructure : public SDependentStructure {
  LayeredHStructure myBase;
  gsl_vector *myInvWeights;
  gsl_vector *myInvSqrtWeights;
  void mulInvWij( gsl_matrix * res, int i ) const;
public:  
  WLayeredHStructure( const double *oldNk, size_t q, int M, 
                     const gsl_vector *weights = NULL );
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
                         bool scaled = true );
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

class WMosaicHStructure : public StripedStructure {
protected:
  static Structure **allocStripe( gsl_vector *oldNk, gsl_vector *oldMl,  
                gsl_vector *Wk );
public:
  WMosaicHStructure( gsl_vector *oldNk, gsl_vector *oldMl, gsl_vector *Wk );
  virtual ~WMosaicHStructure() {}
};





