


class slraLayeredHankelStructure : virtual public slraStructure, virtual public slraWkInterface {
  int myQ;	                /* number of layers */
  size_t myM;
  typedef struct {
    size_t blocks_in_row;       /* Number of blocks in a row of Ci */
    double inv_w;            /* Square root of inverse of the weight */
  } slraFlexBlock;

  
  size_t myNplusD;
  size_t myMaxLag;

  void computeStats();
  void computeWkParams(); 

  gsl_matrix **myA;

  slraFlexBlock *mySA;	/* q-element array describing C1,...,Cq; */  
public:
  slraLayeredHankelStructure( const slraLayeredHankelStructure &s ); /* Copy constructor */
  slraLayeredHankelStructure( const double *oldNk, size_t q, int M, 
                     const double *w_k = NULL );
  virtual ~slraLayeredHankelStructure();

  virtual int getNplusD() const { return myNplusD; }
  virtual int getM() const { return myM; }
  virtual int getNp() const { return (myM - 1) * myQ + myNplusD; }
  virtual slraGammaCholesky *createGammaComputations( int r, double reg_gamma ) const;
  virtual slraDGamma *createDerivativeComputations( int r ) const;

  void setM( int m );
  
  int getQ() const { return myQ; }
 
  int getMaxLag() const { return myMaxLag; }
  
  
  int getFlexBlockLag( int l ) const { return mySA[l].blocks_in_row; }
  int getFlexBlockNCol( int l ) const { return mySA[l].blocks_in_row; }

  bool isFlexBlockExact( int l ) const { return (mySA[l].inv_w == 0.0); }
  double getInvBlockWeight( int l ) const { return mySA[l].inv_w; }
  
  int getFlexBlockNp( int l ) const { return getFlexBlockLag(l) + getM() - 1; }
  
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ); 
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );

  virtual int getS() const { return myMaxLag; }
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





