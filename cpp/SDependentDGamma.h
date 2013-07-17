class SDependentDGamma : public DGamma {
private:
  const SDependentStructure *myW;
  size_t myD;
  gsl_vector *myTmp1, *myTmp2, * myTmp3;
  gsl_vector *myYrR;
  gsl_matrix *myEye;
public:  
  SDependentDGamma( const SDependentStructure *s, size_t D );
  virtual ~SDependentDGamma();
  virtual void calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                   const gsl_vector *yr );
  virtual void calcDijGammaYr( gsl_vector *res, const gsl_matrix *R, 
                   size_t i, size_t j, const gsl_vector *Yr );
};
