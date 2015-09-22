/** The core class for cost function/derivatives computation.
 * This is the class for compuation of the VARPRO cost function, 
 * its gradient and Jacobian, according to \cite slra-efficient.
 *
 * The matrix structure that is allowed is
 */
class VarproFunction  {
  Structure *myStruct;
  size_t myD;
  Cholesky *myGam;
  DGamma *myDeriv;
  double myReggamma;
  bool myIsGCD;

  gsl_matrix *myPhi;
  gsl_matrix *myMatr;
  gsl_matrix *myRorig;
  gsl_matrix *myTmpGradR, *myTmpGradR2;
  gsl_matrix *myTmpJac, *myTmpJac2, *myTmpJtJ, *myTmpEye;
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  gsl_vector *myP;
  double myPWnorm2;

  gsl_vector *myPhiPermCol;  
  gsl_vector *myTmpJacobianCol;  
protected:  
  void setPhiPermCol( size_t i, const gsl_matrix *perm, gsl_vector *phiPermCol );
  virtual void fillZmatTmpJac( gsl_matrix *Zmatr, const gsl_vector* yr,
                               const gsl_matrix *PhiTRt, double factor = 0.5 );
  virtual void mulZmatPerm( gsl_vector* res, const gsl_matrix *Zmatr,
                            const gsl_matrix *perm, size_t i, size_t j );
  size_t getM() { return myStruct->getM(); }
  
  virtual void computeGammaSr( const gsl_matrix *R, gsl_matrix *PhiTRt,
                               gsl_vector *Sr, bool regularize_gamma );
  virtual void computePseudoJacobianLsFromYr( const gsl_vector* yr, 
                   const gsl_matrix *Rorig, const gsl_matrix *perm, 
                   gsl_matrix *pjac, double factor = 0.5 );

  virtual void computeJacobianOfCorrection( const gsl_vector* yr, 
                   const gsl_matrix *Rorig, const gsl_matrix *perm, gsl_matrix *jac );
  virtual void computeGradFromYr( const gsl_vector* yr, const gsl_matrix *Rorig, 
                                  const gsl_matrix *perm, gsl_matrix *grad );
  const gsl_vector *getP() { return myP; }
  virtual const gsl_matrix * getOrigSMatr() { return myMatr; }
public:
  
  
  VarproFunction( const gsl_vector *p, Structure *s,
                  size_t d, gsl_matrix *PhiT, bool isGCD = false );
  virtual ~VarproFunction();
  
  bool isGCD() { return myIsGCD; }

  size_t getD() { return myD; }
  size_t getN() { return myStruct->getN(); }
  size_t getNrow() { return myPhi->size2; }
  size_t getNp() { return myStruct->getNp(); }
  double getReggamma() { return myReggamma; }
  void setReggamma( double reg_gamma ) { myReggamma = reg_gamma; }


  virtual void computeFuncAndGrad( const gsl_matrix* R, double* f, 
                                   const gsl_matrix *perm, gsl_matrix *gradR );
  virtual void computePhat( gsl_vector* p, const gsl_matrix* R );
  virtual void computeCorrectionAndJacobian( const gsl_matrix* R, 
                   const gsl_matrix *perm, gsl_vector *res, gsl_matrix *jac  );

  void computeDefaultRTheta( gsl_matrix *RTheta ); 
  virtual void computeFuncAndPseudoJacobianLs( const gsl_matrix* R, gsl_matrix *perm,
                   gsl_vector *res, gsl_matrix *jac, double factor = 0.5 );
  virtual void computeJtJmulE( const gsl_matrix* R, const gsl_matrix* E,  gsl_matrix *out, int useJtJ = 1 );
};




