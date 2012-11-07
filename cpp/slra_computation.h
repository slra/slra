class CostFunction {
  Structure *myStruct;
  int myD;
  Cholesky *myGam;
  DGamma *myDeriv;
  gsl_matrix *myMatr;
  gsl_matrix *myPerm;
  
  gsl_matrix *myTmpThetaExt;
  gsl_matrix *myTmpGradR;
  gsl_matrix *myTmpGradR2;
  gsl_matrix *myTmpR;  
  gsl_matrix *myEye;  
  gsl_matrix *myTmpJac;  
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  const gsl_vector *myP;
  double myPNorm;

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
  void computeRGammaSr( const gsl_vector *x, gsl_matrix *R, gsl_vector *Sr, 
          bool regularize_gamma = true ) ;

  void computeJacobianZij( gsl_vector*res, int i, int j, gsl_vector* yr, 
                           gsl_matrix *R, double factor = 0.5 );
  void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *R, 
                                      gsl_matrix *jac );

  void computeJacobianOfCorrection( gsl_vector* yr, gsl_matrix *R, 
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
