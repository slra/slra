class SDependentCholesky : public Cholesky {
protected:
  const SDependentStructure *myW;
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
  SDependentCholesky( const SDependentStructure *s, int D,  
                            double reg_gamma );
  virtual ~SDependentCholesky();

  int getD() const { return myD; }
  int getM() const { return myW->getM(); }
  int getNplusD() const { return myW->getNplusD(); }
  int getS() const { return myW->getS(); }


  virtual void calcGammaCholesky( gsl_matrix *R );

  
  virtual void multInvPartCholeskyArray( double * yr, int trans, 
                   size_t size, size_t chol_size );
  virtual void multInvPartGammaArray( double * yr, size_t size, 
                   size_t chol_size );

  virtual void multInvCholeskyVector( gsl_vector * yr, int trans );
  virtual void multInvGammaVector( gsl_vector * yr );
  virtual void multInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans );
};

class StationaryCholesky : public SDependentCholesky {
protected:
  const StationaryStructure *myWs;
  gsl_matrix *myGamma;
  gsl_matrix *myWkTmp;
  
  virtual void computeGammak( gsl_matrix *R );
public:
  StationaryCholesky( const StationaryStructure *s, int D, 
                             double reg_gamma  );
  virtual ~StationaryCholesky();

  virtual void computeGammaUpperPart( gsl_matrix *R );
};


#ifdef USE_SLICOT
class StationaryCholeskySlicot : public StationaryCholesky {
  double *myGammaVec;
  double *myCholeskyWork;
  size_t myCholeskyWorkSize;
public:
  StationaryCholeskySlicot( const StationaryStructure *s, int D, 
                            double reg_gamma );
  virtual ~StationaryCholeskySlicot();

  virtual void computeCholeskyOfGamma( gsl_matrix *R );
};
#endif /* USE_SLICOT */


class SameStripedStationaryCholesky : public Cholesky {
  StationaryCholesky *myBase;
  const MosaicHStructure *myStruct;
public:  
  SameStripedStationaryCholesky( const MosaicHStructure *s, 
      int r, int use_slicot, double reg_gamma );
  virtual ~SameStripedStationaryCholesky();

  virtual void calcGammaCholesky( gsl_matrix *R );
  virtual void multInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multInvGammaVector( gsl_vector * yr );                
};




