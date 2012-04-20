class slraDGammaBTBanded : public slraDGamma {
  const slraStationaryStructure *myW;
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
  slraDGammaBTBanded( const slraStationaryStructure *s, int r );
  virtual ~slraDGammaBTBanded();
  int getD() const { return myD; }

  
  virtual void calcYrtDgammaYr( gsl_matrix *mgrad_r, gsl_matrix *R, gsl_vector *yr );
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );

};


class slraDGammaBBanded : public slraDGamma {
  const slraSDependentStructure *myW;
  size_t myD;

public:  
  slraDGammaBBanded( const slraSDependentStructure *s, int r );
  virtual ~slraDGammaBBanded();
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr );
  virtual void calcYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, 
                   gsl_vector *yr ) {}
};
