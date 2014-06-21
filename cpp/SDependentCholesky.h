class MuDependentCholesky : public Cholesky {
protected:
  const MuDependentStructure *myW;
  size_t myD;
  double my_reg_gamma;
    
  size_t myMu_1;
  size_t myDMu;
  size_t myDN;
  size_t myDMu_1;
  gsl_matrix *myTempR, *myTempVktR, *myTempGammaij;
  
  double *myPackedCholesky;
protected:  
  virtual void computeGammaUpperPart( const gsl_matrix *R, double reg = 0 );

public:
  MuDependentCholesky( const MuDependentStructure *s, size_t D );
  virtual ~MuDependentCholesky();

  size_t getD() const { return myD; }
  size_t getN() const { return myW->getN(); }
  size_t getM() const { return myW->getM(); }
  size_t getMu() const { return myW->getMu(); }

  virtual void calcGammaCholesky( const gsl_matrix *R, double reg = 0 );

  virtual void multInvCholeskyVector( gsl_vector * yr, long trans );
  virtual void multInvGammaVector( gsl_vector * yr );
};








