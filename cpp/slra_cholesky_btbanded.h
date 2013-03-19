class SDependentCholesky : public Cholesky {
protected:
  const SDependentStructure *myW;
  size_t myD;
  double my_reg_gamma;
    
  size_t myS_1;
  size_t myDS;
  size_t myDN;
  size_t myDS_1;
  gsl_matrix *myTempR, *myTempWktR, *myTempGammaij;
  
  double *myPackedCholesky;
protected:  
  virtual void computeGammaUpperPart( const gsl_matrix *R, double reg = 0 );

public:
  SDependentCholesky( const SDependentStructure *s, size_t D );
  virtual ~SDependentCholesky();

  size_t getD() const { return myD; }
  size_t getN() const { return myW->getN(); }
  size_t getM() const { return myW->getM(); }
  size_t getS() const { return myW->getS(); }

  virtual void calcGammaCholesky( const gsl_matrix *R, double reg = 0 );

  virtual void multInvCholeskyVector( gsl_vector * yr, long trans );
  virtual void multInvGammaVector( gsl_vector * yr );
};








