class slraGammaCholeskyBBanded : public slraGammaCholesky {
protected:
  const slraStructure *myW;
  size_t myD;
  double my_reg_gamma;
    
  size_t s_minus_1;
  size_t d_times_s;
  size_t d_times_Mg;
  size_t d_times_s_minus_1;
  
  double *myPackedCholesky;
public:
  slraGammaCholeskyBBanded( const slraStructure *s, int r, double reg_gamma );
  virtual ~slraGammaCholeskyBBanded();

  int getD() const { return myD; }
  int getM() const { return myW->getM(); }
  int getNplusD() const { return myW->getNplusD(); }
  int getS() const { return myW->getS(); }
  
  virtual void multiplyInvPartCholeskyArray( double * yr, int trans, size_t size, size_t chol_size );
  virtual void multiplyInvPartGammaArray( double * yr, size_t size, size_t chol_size );

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans );
  virtual void multiplyInvGammaVector( gsl_vector * yr );
  virtual void multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans );
};

class slraGammaCholeskyBTBanded : public slraGammaCholeskyBBanded {
  const slraStationaryStructure *myWs;

  gsl_matrix *myGamma;
  gsl_matrix *myWkTmp;
#ifdef USE_SLICOT
  double *myGammaVec;
  double *myCholeskyWork;
  int myCholeskyWorkSize;
#endif  
public:
  slraGammaCholeskyBTBanded( const slraStationaryStructure *s, int r, double reg_gamma  );
  virtual ~slraGammaCholeskyBTBanded();

  virtual void computeCholeskyOfGamma( gsl_matrix *R );
};

class slraGammaCholeskySameDiagBTBanded : public slraGammaCholesky {
  slraGammaCholeskyBTBanded *myBase;
  const slraMosaicHankelStructure *myStruct;
public:  
  slraGammaCholeskySameDiagBTBanded( const slraMosaicHankelStructure *s, int r, int use_slicot, double reg_gamma  ) :
      myStruct(s) {
    myBase = (slraGammaCholeskyBTBanded *)myStruct->getMaxBlock()->createGammaComputations(r, reg_gamma);  
  }
  virtual ~slraGammaCholeskySameDiagBTBanded() {
    delete myBase;
  }
  
  virtual void computeCholeskyOfGamma( gsl_matrix *R ) {
    myBase->computeCholeskyOfGamma(R);
  }
  
  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multiplyInvGammaVector( gsl_vector * yr );                
};



class slraGammaCholeskyBBandedLH : public slraGammaCholeskyBBanded {
  const slraLayeredHankelWeightedStructure *myWs;
  gsl_matrix *myTempR, *myTempRtWkt, *myTempGammaij;
  gsl_vector *myTempInvW;
public:
  slraGammaCholeskyBBandedLH( const slraLayeredHankelWeightedStructure *s, int r, double reg_gamma ) :
      slraGammaCholeskyBBanded(s, r, reg_gamma), myWs(s) {
    myTempR = gsl_matrix_alloc(getNplusD(), getD());
    myTempRtWkt = gsl_matrix_alloc(getD(), getNplusD());
    myTempGammaij = gsl_matrix_alloc(getD(), getD());
    myTempInvW = gsl_vector_alloc(getNplusD());
  }
  virtual ~slraGammaCholeskyBBandedLH() {
    gsl_matrix_free(myTempR);
    gsl_matrix_free(myTempRtWkt);
    gsl_matrix_free(myTempGammaij);
    gsl_vector_free(myTempInvW);
  }

  virtual void computeCholeskyOfGamma( gsl_matrix *R );
};
