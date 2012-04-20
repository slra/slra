class StripedStructure : public slraStructure {
  size_t myN;

  /* Helper variables */
  size_t myM;
  size_t myNp;

  size_t myMaxMlInd;
protected:
  StripedStructure( size_t N, slraStructure **stripe );

public:
  slraStructure **myStripe;
  virtual ~StripedStructure();

  virtual int getM() const { return myM; }
  virtual int getNp() const { return myNp; }
  virtual int getNplusD() const { return myStripe[0]->getNplusD(); }
  
  int getBlocksN() const { return myN; }

  int getMl( int k ) const { return myStripe[k]->getM(); }
  const slraStructure *getBlock( size_t k ) const { 
    return myStripe[k]; 
  }

  int getMaxMl() const { return getMl(myMaxMlInd); }
  const slraStructure *getMaxBlock() const { 
    return getBlock(myMaxMlInd); 
  }

  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) ;
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );

  virtual slraGammaCholesky *createGammaComputations( int r, double reg_gamma ) const;
  virtual slraDGamma *createDerivativeComputations( int r ) const;
};

class StripedCholesky : virtual public slraGammaCholesky {
  slraGammaCholesky **myGamma;
  int myD;
  const StripedStructure *myStruct;
public:  
  StripedCholesky( const StripedStructure *s, int r, double reg_gamma );
  virtual ~StripedCholesky();

  
  virtual void computeCholeskyOfGamma( gsl_matrix *R );

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multiplyInvGammaVector( gsl_vector * yr );                
};

class StripedDGamma : virtual public slraDGamma {
  slraDGamma **myLHDGamma;
  const StripedStructure *myStruct;
  gsl_matrix *myTmpGrad;
public:  
  StripedDGamma( const StripedStructure *s, int r  ) ;
  virtual ~StripedDGamma();

  virtual void calcYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, gsl_vector *yr );

  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );

};

