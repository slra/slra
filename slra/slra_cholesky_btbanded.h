class slraGammaCholeskyBBanded : public slraGammaCholesky {
protected:
  const slraSDependentStructure *myW;
  size_t myD;
  double my_reg_gamma;
    
  size_t s_minus_1;
  size_t d_times_s;
  size_t d_times_Mg;
  size_t d_times_s_minus_1;
  gsl_matrix *myTempR, *myTempWktR, *myTempGammaij;
  
  double *myPackedCholesky;

protected:  
  virtual void computeGammaUpperPart( gsl_matrix *R );

public:
  slraGammaCholeskyBBanded( const slraSDependentStructure *s, int r, double reg_gamma );
  virtual ~slraGammaCholeskyBBanded();

  int getD() const { return myD; }
  int getM() const { return myW->getM(); }
  int getNplusD() const { return myW->getNplusD(); }
  int getS() const { return myW->getS(); }


  virtual void computeCholeskyOfGamma( gsl_matrix *R );

  
  virtual void multiplyInvPartCholeskyArray( double * yr, int trans, size_t size, size_t chol_size );
  virtual void multiplyInvPartGammaArray( double * yr, size_t size, size_t chol_size );

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans );
  virtual void multiplyInvGammaVector( gsl_vector * yr );
  virtual void multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans );
};

class slraGammaCholeskyBTBanded : public slraGammaCholeskyBBanded {
protected:
  const slraStationaryStructure *myWs;

  gsl_matrix *myGamma;
  gsl_matrix *myWkTmp;
  
  virtual void computeGammak( gsl_matrix *R );
public:
  slraGammaCholeskyBTBanded( const slraStationaryStructure *s, int r, double reg_gamma  );
  virtual ~slraGammaCholeskyBTBanded();

  virtual void computeGammaUpperPart( gsl_matrix *R );
};


#ifdef USE_SLICOT
class slraGammaCholeskyBTBandedSlicot : public slraGammaCholeskyBTBanded {
  double *myGammaVec;
  double *myCholeskyWork;
  int myCholeskyWorkSize;
public:
  slraGammaCholeskyBTBanded( const slraStationaryStructure *s, int r, double reg_gamma  );
  virtual ~slraGammaCholeskyBTBanded();

  virtual void computeGammaUpperPart( gsl_matrix *R );
};
#endif /* USE_SLICOT */


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




