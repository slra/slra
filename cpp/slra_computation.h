
class CostFunction  {
  Structure *myStruct;
  size_t myD;
  Cholesky *myGam;
  DGamma *myDeriv;
  double myReggamma;

  gsl_matrix *myPhi;
  gsl_matrix *myMatr;
  gsl_matrix *myRorig;
  gsl_matrix *myTmpGradR, *myTmpGradR2;
  gsl_matrix *myTmpJac, *myTmpJac2;  
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  gsl_vector *myP;

  gsl_vector *myPhiPermCol;  
  gsl_vector *myTmpJacobianCol;  
protected:  
  void setPhiPermCol( size_t i, gsl_matrix *perm );
  virtual void computeZmatTmpJac( gsl_vector* yr, gsl_matrix *Rorig, double factor = 0.5 );
  virtual void mulZmatPerm( gsl_vector* res, gsl_matrix *perm, size_t i, size_t j );

  size_t getM() { return myStruct->getM(); }
  
  virtual void computeGammaSr( const gsl_matrix *R, gsl_matrix *Rorig, gsl_vector *Sr, 
                               bool regularize_gamma );
  virtual void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *Rorig, 
                       gsl_matrix *perm, gsl_matrix *jac, double factor = 0.5 );
  virtual void computeJacobianOfCorrection( gsl_vector* yr, gsl_matrix *Rorig, 
                                    gsl_matrix *perm, gsl_matrix *jac );
  virtual void computeGradFromYr( gsl_vector* yr, const gsl_matrix *Rorig, 
                                  gsl_matrix *perm, gsl_matrix *grad );
public:
  CostFunction( const gsl_vector *p, Structure *s, size_t d, gsl_matrix *Phi );
  virtual ~CostFunction();
  
  size_t getD() { return myD; }
  size_t getN() { return myStruct->getN(); }
  size_t getNrow() { return myPhi->size2; }
  size_t getNp() { return myStruct->getNp(); }
  const gsl_vector *getP() { return myP; }
  double getReggamma() { return myReggamma; }
  void setReggamma( double reg_gamma ) { myReggamma = reg_gamma; }

  void computeDefaultRTheta( gsl_matrix *RTheta ); 

  virtual void computePhat( gsl_vector* p, gsl_matrix* R );
  virtual void computeCorrectionAndJacobian( gsl_matrix* R, gsl_matrix *perm,
                                    gsl_vector *res, gsl_matrix *jac  );
  virtual void computeFuncAndPseudoJacobianLs( const gsl_matrix* R, gsl_matrix *perm,
                   gsl_vector *res, gsl_matrix *jac, double factor = 0.5 );
  virtual void computeFuncAndGrad( const gsl_matrix* R, double* f, 
                                   gsl_matrix *perm, gsl_matrix *gradR );
  virtual void multiplyByHessian( const gsl_matrix* R, const gsl_matrix* Hin, 
                                  gsl_matrix* Hout );
  virtual void computeHessian( const gsl_matrix* R, gsl_matrix* H );

  virtual const gsl_matrix * getOrigSMatr() { return myMatr; }
};




