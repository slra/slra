class StationaryCholesky : public SDependentCholesky {
public:
  StationaryCholesky( const StationaryStructure *s, size_t D );
  virtual ~StationaryCholesky();

  virtual void computeGammaUpperPart( const gsl_matrix *R, double reg = 0 );

protected:
  const StationaryStructure *myWs;
  gsl_matrix *myGamma;
  gsl_matrix *myWkTmp;
  
  virtual void computeGammak( const gsl_matrix *R, double reg = 0 );
};


