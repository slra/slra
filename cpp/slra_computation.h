
class CostFunction  {
  Structure *myStruct;
  int myD;
  Cholesky *myGam;
  DGamma *myDeriv;

  gsl_matrix *myPhi;
  gsl_matrix *myMatr;
  gsl_matrix *myRorig;
  gsl_matrix *myTmpGradR, *myTmpGradR2;
  gsl_matrix *myTmpJac;  
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  const gsl_vector *myP;

  gsl_vector *myPhiPermCol;  
  gsl_vector *myTmpJacobianCol;  
  gsl_matrix *myTmpGrad;  
protected:  
  virtual void computeZmatTmpJac( gsl_vector* yr, gsl_matrix *Rorig, double factor = 0.5 );
  virtual void mulZmatPerm( gsl_vector* res, gsl_matrix *perm, size_t i, size_t j );

  int getM() { return myStruct->getM(); }
  
  virtual void computeGammaSr( const gsl_matrix *R, gsl_matrix *Rorig, gsl_vector *Sr, 
                               bool regularize_gamma );
  virtual void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *Rorig, 
                                      gsl_matrix *perm, gsl_matrix *jac );
  virtual void computeJacobianOfCorrection( gsl_vector* yr, gsl_matrix *Rorig, 
                                    gsl_matrix *perm, gsl_matrix *jac );
  virtual void computeGradFromYr( gsl_vector* yr, const gsl_matrix *Rorig, 
                                  gsl_matrix *perm, gsl_matrix *grad );
public:
  CostFunction( Structure *s, int d, const gsl_vector *p, double reggamma, gsl_matrix *Phi );
  virtual ~CostFunction();
  
  int getD() { return myD; }
  int getN() { return myStruct->getN(); }
  int getNrow() { return myPhi->size2; }
  int getNp() { return myStruct->getNp(); }

  void computeDefaultRTheta( gsl_matrix *RTheta ); 

  virtual void computeCorrection( gsl_vector* p, gsl_matrix* R );
  virtual void computeCorrectionAndJacobian( gsl_matrix* R, gsl_matrix *perm,
                                    gsl_vector *res, gsl_matrix *jac  );
  virtual void computeFuncAndPseudoJacobianLs( gsl_matrix* R, gsl_matrix *perm,
                                       gsl_vector *res, gsl_matrix *jac );
  virtual void computeFuncAndGrad( const gsl_matrix* R, double* f, 
                                   gsl_matrix *perm, gsl_matrix *gradR );

  virtual const gsl_matrix * getOrigSMatr() { return myMatr; }
};




