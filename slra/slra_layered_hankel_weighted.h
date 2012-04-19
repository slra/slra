class slraLayeredHankelWeightedStructure : public slraSDependentStructure {
  slraLayeredHankelStructure myBase;
  gsl_vector *myInvSqrtWeights;
public:
  slraLayeredHankelWeightedStructure( const double *oldNk, size_t q, int M, 
                     const double *weights = NULL, bool block_weights = false );
  virtual ~slraLayeredHankelWeightedStructure();

  /* slraLayeredHankelStructure delegated methods */
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

  /* slraStructure methods */
  virtual slraGammaCholesky *createGammaComputations( int r, double reg_gamma ) const;
  virtual slraDGamma *createDerivativeComputations( int r ) const;
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );

  /* slraSDependentStructure methods */
  virtual void WijB( gsl_matrix *res, int i, int j, const gsl_matrix *B ) const;
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const;
  virtual void AtWijV( gsl_vector *res, int i, int j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const;
                      
  /* Structure-specific methods */
  double getInvSqrtWeights( int i ) const { return gsl_vector_get(myInvSqrtWeights, i); }
  void mulInvWij( gsl_matrix * res, int i ) const;
};




