

class slraLayeredHankelStructure : public slraStationaryStructure {
  int myQ;	                /* number of layers */
  size_t myM;
  typedef struct {
    size_t blocks_in_row;       /* Number of blocks in a row of Ci */
    double inv_w;            /* Square root of inverse of the weight */
  } slraLayer;

  
  size_t myNplusD;
  size_t myMaxLag;

  void computeStats();
  void computeWkParams(); 

  gsl_matrix **myA;

  slraLayer *mySA;	/* q-element array describing C1,...,Cq; */  
protected:
  int nvGetNp() const { return (myM - 1) * myQ + myNplusD; }  
  
public:
  slraLayeredHankelStructure( const double *oldNk, size_t q, int M, 
                     const double *layer_w = NULL );
  virtual ~slraLayeredHankelStructure();

  virtual int getNplusD() const { return myNplusD; }
  virtual int getM() const { return myM; }
  virtual int getS() const { return myMaxLag; }
  
  virtual int getNp() const { return nvGetNp(); }
  virtual slraGammaCholesky *createGammaComputations( int r, double reg_gamma ) const;
  virtual slraDGamma *createDerivativeComputations( int r ) const;

  void setM( int m );
  
  int getQ() const { return myQ; }
  int getMaxLag() const { return myMaxLag; }
  
  int getLayerLag( int l ) const { return mySA[l].blocks_in_row; }
  bool isLayerExact( int l ) const { return (mySA[l].inv_w == 0.0); }
  double getLayerInvWeight( int l ) const { return mySA[l].inv_w; }
  int getLayerNp( int l ) const { return getLayerLag(l) + getM() - 1; }
  
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ); 
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );


  virtual const gsl_matrix *getWk( int k ) const { 
    return myA[k]; 
  }
};


class slraMosaicHankelStructure : public slraStripedStructure {
public:
  static slraStructure **allocStripe( size_t q, size_t N, double *oldNk, double *oldMl, double *Wk );

  slraMosaicHankelStructure( size_t q, size_t N, double *oldNk, double *oldMl, double *Wk ) :
      slraStripedStructure(N, allocStripe(q, N, oldNk, oldMl, Wk)) {
  }
  virtual ~slraMosaicHankelStructure() {}

  virtual slraGammaCholesky *createGammaComputations( int r, double reg_gamma ) const;
};


class slraLayeredHankelWeightedStructure : public slraLayeredHankelStructure {
  gsl_vector *myInvSqrtWeights;
public:
  slraLayeredHankelWeightedStructure( const double *oldNk, size_t q, int M, 
                     const double *weights = NULL );
  virtual ~slraLayeredHankelWeightedStructure();

  double getInvSqrtWeights( int i ) const { return gsl_vector_get(myInvSqrtWeights, i); }

  virtual slraGammaCholesky *createGammaComputations( int r, double reg_gamma ) const;
  virtual slraDGamma *createDerivativeComputations( int r ) const;
  
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );
};




