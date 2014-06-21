class StationaryDGamma : public DGamma {
private:
  const StationaryStructure *myW;
  size_t  myD;
  
  gsl_vector *myTempVkColRow;
  gsl_vector *myDGammaVec;
  gsl_matrix *myDGammaTrMat;
  gsl_matrix *myDGamma;
  gsl_vector *myTmpCol;
  
  gsl_matrix *myVk_R;
  gsl_matrix *myN_k;
  gsl_matrix *myEye;
public:
  StationaryDGamma( const StationaryStructure *s, size_t D );
  virtual ~StationaryDGamma();
  size_t getD() const { return myD; }

  
  virtual void calcYrtDgammaYr( gsl_matrix *mgrad_r, const gsl_matrix *R, 
                   const gsl_vector *yr );
  virtual void calcDijGammaYr( gsl_vector *res, const gsl_matrix *R, 
                    size_t i, size_t j, const gsl_vector *Yr );
};



