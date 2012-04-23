class StationaryDGamma : public DGamma {
  const StationaryStructure *myW;
  size_t  myD;
  
  gsl_vector *myTempWkColRow;
  gsl_vector *myDGammaVec;
  gsl_matrix *myDGammaTrMat;
  gsl_matrix *myDGamma;
  gsl_vector *myTmpCol;
  
  gsl_matrix *myWk_R;
  gsl_matrix *myWkT_R;
  gsl_matrix *myN_k;

public:
  StationaryDGamma( const StationaryStructure *s, int D );
  virtual ~StationaryDGamma();
  int getD() const { return myD; }

  
  virtual void calcYrtDgammaYr( gsl_matrix *mgrad_r, gsl_matrix *R, gsl_vector *yr );
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );

};


class SDependentDGamma : public DGamma {
  const SDependentStructure *myW;
  size_t myD;
  gsl_vector *myTmp1, *myTmp2;
  gsl_vector *myYrR;
public:  
  SDependentDGamma( const SDependentStructure *s, int D );
  virtual ~SDependentDGamma();
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );
  virtual void calcYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, 
                   gsl_vector *yr ) {}
};
