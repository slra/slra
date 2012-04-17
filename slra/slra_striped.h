class slraStripedStructure : public slraStructure {
  slraStructure **myLHStripe;
  size_t myN;

  /* Helper variables */
  size_t myM;
  size_t myNp;

  size_t myMaxMlInd;
protected:
  slraStripedStructure( size_t N, slraStructure **stripe );

public:
  virtual ~slraStripedStructure();

  virtual int getM() const { return myM; }
  virtual int getNp() const { return myNp; }
  virtual int getNplusD() const { return myLHStripe[0]->getNplusD(); }
  virtual int getS() const { return myLHStripe[0]->getS(); }
  
  int getBlocksN() const { return myN; }

  int getMl( int k ) const { return myLHStripe[k]->getM(); }
  const slraStructure *getBlock( size_t k ) const { 
    return myLHStripe[k]; 
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

class slraDiagGammaCholesky : virtual public slraGammaCholesky {
  slraGammaCholesky **myGamma;
  int myD;
  const slraStripedStructure *myStruct;
public:  
  slraDiagGammaCholesky( const slraStripedStructure *s, int r, double reg_gamma  ) ;
  virtual ~slraDiagGammaCholesky();

  
  virtual void computeCholeskyOfGamma( gsl_matrix *R );

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multiplyInvGammaVector( gsl_vector * yr );                
};

class slraDGammaStriped : virtual public slraDGamma {
  slraDGamma **myLHDGamma;
  const slraStripedStructure *myStruct;
  gsl_matrix *myTmpGrad;
public:  
  slraDGammaStriped( const slraStripedStructure *s, int r  ) ;
  virtual ~slraDGammaStriped();

  virtual void computeYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, gsl_vector *yr );

  virtual void computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );

};

