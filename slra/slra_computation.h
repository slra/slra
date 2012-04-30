class CostFunction {
  Structure *myStruct;
  int myRank;
  Cholesky *myGam;
  DGamma *myDeriv;
  gsl_matrix *myMatr;
  gsl_matrix *myPerm;
  gsl_matrix *myMatrMulPerm;
  
  gsl_matrix *myTmpThetaExt;
  gsl_matrix *myTmpGradR;
  gsl_matrix *myTmpGradR2;
  gsl_matrix *myTmpR;  
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  const gsl_vector *myP;
  double myPNorm;
  bool isGCD;

  gsl_vector *myTmpJacobianCol;  
  gsl_matrix *myTmpGrad;  
public:
  CostFunction( Structure *s, int r, const gsl_vector *p, 
                    opt_and_info *opt, gsl_matrix *perm  );
  virtual ~CostFunction();
  
  int getD() { return myPerm->size2 - myRank; }
  int getNplusD() { return myStruct->getNplusD(); }
  int getRank() { return myRank; }
  int getRsize() { return myPerm->size2; }
  int getM() { return myStruct->getM(); }
  int getNp() { return myStruct->getNp(); }
  const gsl_matrix * getPerm() { return myPerm; }
  const gsl_matrix * getSMart() { return myMatr; }
  const gsl_matrix * getPhiSMatr() { return myMatrMulPerm; }

  void computeRTheta( const gsl_matrix * x_mat, gsl_matrix *RTheta ); 
  void computeDefaultRTheta( gsl_matrix *RTheta ); 
  static void R_2_x( const gsl_matrix *R, gsl_matrix * x );
  
  void computeR( gsl_matrix_const_view x_mat, gsl_matrix *R ); 
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
  void computeSr( gsl_matrix *R, gsl_vector *Sr );
  void computeRGammaSr( const gsl_vector *x, gsl_matrix *R, gsl_vector *Sr ) ;

  void computeJacobianZij( gsl_vector*res, int i, int j, gsl_vector* yr, 
                           gsl_matrix *R, double factor = 0.5 );
  void computePseudoJacobianCorrectFromYr( gsl_vector* yr, gsl_matrix *R, 
                                           gsl_matrix *jac );
  void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *R, 
                                      gsl_matrix *jac );
  void computeGradFromYr( gsl_vector* yr, gsl_matrix *R, gsl_vector *grad );
  void computeCorrectionAndJacobian( const gsl_vector* x, gsl_vector *res, 
                                     gsl_matrix *jac  );
  void computeFuncAndPseudoJacobianLs( const gsl_vector* x, gsl_vector *res, 
                                       gsl_matrix *jac );
  void computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad );
  
  void  computeCorrection( gsl_vector* p, const gsl_vector* x );
  
  static int _f_ls( const gsl_vector* x, void* params, gsl_vector* res ) {
    ((CostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int _df_ls( const gsl_vector* x,  void* params, gsl_matrix* jac ) {
    ((CostFunction *)params)->computeFuncAndPseudoJacobianLs(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int _fdf_ls( const gsl_vector* x,  void* params, gsl_vector* res,
                      gsl_matrix *jac ) {
    ((CostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, jac);
    return GSL_SUCCESS;
  }

  static int _f_cor( const gsl_vector* x, void* params, gsl_vector* res ) {
    ((CostFunction *)params)->computeCorrectionAndJacobian(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int _df_cor( const gsl_vector* x,  void* params, gsl_matrix* jac ) {
    ((CostFunction *)params)->computeCorrectionAndJacobian(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int _fdf_cor( const gsl_vector* x, void* params, gsl_vector*res, 
                       gsl_matrix *jac ) {
    ((CostFunction *)params)->computeCorrectionAndJacobian(x, res, jac);
    return GSL_SUCCESS;
  }

  static double _f( const gsl_vector* x, void* params ) {
    double f;
    ((CostFunction *)params)->computeFuncAndGrad(x, &f, NULL);
    return f;
  }
  static void _df( const gsl_vector* x, void* params, gsl_vector *grad ) {
    ((CostFunction *)params)->computeFuncAndGrad(x, NULL, grad);
  }
  static void _fdf( const gsl_vector* x, void* params, double* f, 
                     gsl_vector* grad ) {
    ((CostFunction *)params)->computeFuncAndGrad(x, f, grad);
  }
};
