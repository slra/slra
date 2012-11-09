class CostFunction {
  Structure *myStruct;
  int myD;
  Cholesky *myGam;
  DGamma *myDeriv;
  gsl_matrix *myMatr;
  gsl_matrix *myPerm;
  
  gsl_matrix *myTmpThetaExt;
  gsl_matrix *myTmpGradR;
  gsl_matrix *myTmpR;  
  gsl_matrix *myEye;  
  gsl_matrix *myTmpJac;  
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  const gsl_vector *myP;

  gsl_vector *myTmpJacobianCol;  
  gsl_matrix *myTmpGrad;  
public:
  CostFunction( Structure *s, int d, const gsl_vector *p, gsl_matrix *perm, 
                double reggamma );
  virtual ~CostFunction();
  
  int getD() { return myD; }
  int getM() { return myStruct->getM(); }
  int getRank() { return myPerm->size2 - myD; }
  int getRsize() { return myPerm->size2; }
  int getN() { return myStruct->getN(); }
  int getNp() { return myStruct->getNp(); }
  const gsl_matrix * getPerm() { return myPerm; }
  const gsl_matrix * getSMatr() { return myMatr; }

  void X2Rtheta( const gsl_matrix * x_mat, gsl_matrix *RTheta ); 
  void computeDefaultRTheta( gsl_matrix *RTheta ); 
  static void Rtheta2X( const gsl_matrix *R, gsl_matrix * x );
  
  void computeR( gsl_matrix_const_view x_mat, gsl_matrix *R ); 
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
  void computeSr( gsl_matrix *R, gsl_vector *Sr );
  void computeGammaSr( const gsl_matrix *R, gsl_vector *Sr, bool regularize_gamma = true ) ;

  void computeZmatTmpJac( gsl_vector* yr, gsl_matrix *R, double factor = 0.5 );
  void mulZmatPerm( gsl_vector* res, gsl_matrix *perm, size_t i, size_t j );
                           
  void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *R, 
                                      gsl_matrix *jac );

  void computeJacobianOfCorrection( gsl_vector* yr, gsl_matrix *R, 
                                    gsl_matrix *jac );
  void computeCorrectionAndJacobian( const gsl_vector* x, gsl_vector *res, 
                                     gsl_matrix *jac  );
  void computeFuncAndPseudoJacobianLs( const gsl_vector* x, gsl_vector *res, 
                                       gsl_matrix *jac );

  void computeFuncAndGrad( const gsl_matrix* R, double* f, gsl_matrix *gradR );
  
  void computeCorrection( gsl_vector* p, const gsl_vector* x );
};

class OptFunctionSLRA : public OptFunction {
protected:
  CostFunction &myFun;
  gsl_matrix *myTmpR;  
  gsl_matrix *myTmpGradR;
  gsl_matrix *myPerm;
  gsl_matrix *myTmpThetaExt;
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
  void X2Rtheta( const gsl_matrix *x, gsl_matrix *RTheta );
public: 
  OptFunctionSLRA( CostFunction &fun, gsl_matrix *perm );
  virtual ~OptFunctionSLRA();
  virtual int getNvar(); 
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad );
};


class OptFunctionSLRACholesky : public OptFunctionSLRA {
public: 
  OptFunctionSLRACholesky( CostFunction &fun, gsl_matrix *perm  ) : OptFunctionSLRA(fun, perm) {}
  virtual int getNsq() { return myFun.getN() * myFun.getD(); }
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac ) {
    myFun.computeFuncAndPseudoJacobianLs(x, res, jac); 
  }   
};

class OptFunctionSLRACorrection : public OptFunctionSLRA {
public: 
  OptFunctionSLRACorrection( CostFunction &fun, gsl_matrix *perm ) : OptFunctionSLRA(fun, perm)  {}
  virtual int getNsq() { return myFun.getNp(); }
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac ) {
    myFun.computeCorrectionAndJacobian(x, res, jac); 
  }   
};




