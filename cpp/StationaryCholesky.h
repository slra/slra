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

#ifdef USE_SLICOT
class StationaryCholeskySlicot : public StationaryCholesky {
public:
  StationaryCholeskySlicot( const StationaryStructure *s, size_t D );
  virtual ~StationaryCholeskySlicot();

  virtual void calcGammaCholesky( const gsl_matrix *R, double reg = 0 );

private:
  double *myGammaVec;
  double *myCholeskyWork;
  size_t myCholeskyWorkSize;
};
#endif /* USE_SLICOT */

