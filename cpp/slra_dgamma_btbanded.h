class StationaryDGamma : public DGamma {
private:
  const StationaryStructure *myW;
  size_t  myD;
  
  gsl_vector *myTempWkColRow;
  gsl_vector *myDGammaVec;
  gsl_matrix *myDGammaTrMat;
  gsl_matrix *myDGamma;
  gsl_vector *myTmpCol;
  
  gsl_matrix *myWk_R;
  gsl_matrix *myN_k;
  gsl_matrix *myEye;
public:
  StationaryDGamma( const StationaryStructure *s, size_t D );
  virtual ~StationaryDGamma();
  size_t getD() const { return myD; }

  
  virtual void calcYrtDgammaYr( gsl_matrix *mgrad_r, const gsl_matrix *R, 
                   const gsl_vector *yr );
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                    size_t i, size_t j, gsl_vector *Yr );
};


class SDependentDGamma : public DGamma {
public:  
  SDependentDGamma( const SDependentStructure *s, size_t D );
  virtual ~SDependentDGamma();
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   size_t i, size_t j, gsl_vector *Yr );
  virtual void calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                   const gsl_vector *yr );
                   

private:
  const SDependentStructure *myW;
  size_t myD;
  gsl_vector *myTmp1, *myTmp2, * myTmp3;
  gsl_vector *myYrR;
  gsl_matrix *myEye;
};
