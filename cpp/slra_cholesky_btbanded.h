class SDependentCholesky : public Cholesky {
protected:
  const SDependentStructure *myW;
  size_t myD;
  double my_reg_gamma;
    
  size_t s_minus_1;
  size_t d_times_s;
  size_t d_times_n;
  size_t d_times_s_minus_1;
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

  
  virtual void multInvPartCholeskyArray( double * yr, int trans, 
                   size_t size, size_t chol_size );
  virtual void multInvPartGammaArray( double * yr, size_t size, 
                   size_t chol_size );

  virtual void multInvCholeskyVector( gsl_vector * yr, int trans );
  virtual void multInvGammaVector( gsl_vector * yr );
  virtual void multInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans );
};

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


class SameStripedStationaryCholesky : public Cholesky {
public:  
  SameStripedStationaryCholesky( const MosaicHStructure *s, 
      size_t r, int use_slicot );
  virtual ~SameStripedStationaryCholesky();

  virtual void calcGammaCholesky( const gsl_matrix *R, double reg = 0 );
  virtual void multInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multInvGammaVector( gsl_vector * yr );                
  
private:
  StationaryCholesky *myBase;
  const MosaicHStructure *myS;
};




