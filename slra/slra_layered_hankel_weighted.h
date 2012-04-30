class WLayeredHStructure : public SDependentStructure {
  LayeredHStructure myBase;
  gsl_vector *myInvWeights;
public:  
  WLayeredHStructure( const double *oldNk, size_t q, int M, 
                     const gsl_vector *weights = NULL );
  virtual ~WLayeredHStructure();

  /* LayeredHankelStructure delegated methods */
  virtual int getNplusD() const { return myBase.getNplusD(); }
  virtual int getM() const { return myBase.getM(); }
  virtual int getNp() const { return myBase.getNp(); }
  virtual int getS() const { return myBase.getS(); }
  int getQ() const { return myBase.getQ(); }
  int getMaxLag() const { return myBase.getMaxLag(); }
  int getLayerLag( int l ) const { return myBase.getLayerLag(l); }
  int getLayerNp( int l ) const { return myBase.getLayerNp(l); }

  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) {
    myBase.fillMatrixFromP(c, p); 
  }

  /* Structure methods */
  virtual Cholesky *createCholesky( int D, double reg_gamma ) const;
  virtual DGamma *createDGamma( int D ) const;
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );

  /* SDependentStructure methods */
  virtual void WijB( gsl_matrix *res, int i, int j, 
                     const gsl_matrix *B ) const;
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const;
  virtual void AtWijV( gsl_vector *res, int i, int j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const;
                      
  /* Structure-specific methods */
  double getInvWeights( int i ) const { 
    return gsl_vector_get(myInvWeights, i); 
  }
  void mulInvWij( gsl_matrix * res, int i ) const;
};

class WMosaicHStructure : public StripedStructure {
protected:
  static Structure **allocStripe( gsl_vector *oldNk, gsl_vector *oldMl,  
                gsl_vector *Wk );
public:
  WMosaicHStructure( gsl_vector *oldNk, gsl_vector *oldMl, gsl_vector *Wk );
  virtual ~WMosaicHStructure() {}
};





