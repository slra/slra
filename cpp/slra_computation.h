
class CostFunction {
  Structure *myStruct;
  int myD;
  Cholesky *myGam;
  DGamma *myDeriv;
  gsl_matrix *myMatr;
  
  gsl_matrix *myTmpGradR;
  gsl_matrix *myEye;  
  gsl_matrix *myTmpJac;  
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  const gsl_vector *myP;

  gsl_vector *myTmpJacobianCol;  
  gsl_matrix *myTmpGrad;  
protected:  
  void computeZmatTmpJac( gsl_vector* yr, gsl_matrix *R, double factor = 0.5 );
  void mulZmatPerm( gsl_vector* res, gsl_matrix *perm, size_t i, size_t j );
  
public:
  CostFunction( Structure *s, int d, const gsl_vector *p, double reggamma );
  virtual ~CostFunction();
  
  int getD() { return myD; }
  const Structure * getSt() { return myStruct; }
  const gsl_matrix * getSMatr() { return myMatr; }
  
  void computeGammaSr( const gsl_matrix *R, gsl_vector *Sr, bool regularize_gamma = true ) ;

                           
  void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *R, 
                                      gsl_matrix *perm, gsl_matrix *jac );

  void computeJacobianOfCorrection( gsl_vector* yr, gsl_matrix *R, 
                                    gsl_matrix *perm, gsl_matrix *jac );


  void computeCorrectionAndJacobian( gsl_matrix* R, gsl_matrix *perm,
                                    gsl_vector *res, gsl_matrix *jac  );
  void computeFuncAndPseudoJacobianLs( gsl_matrix* R, gsl_matrix *perm,
                                       gsl_vector *res, gsl_matrix *jac );
  void computeFuncAndGrad( const gsl_matrix* R, double* f, gsl_matrix *perm, gsl_matrix *gradR );
  
  void computeCorrection( gsl_vector* p, gsl_matrix* R );
};





