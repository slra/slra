class slraCostFunction {
  slraStructure *myStruct;
  int myRank;
  slraGammaCholesky *myGam;
  slraDGamma *myDeriv;
  gsl_matrix *myMatr;
  gsl_matrix *myPerm;
  
  gsl_matrix *myMatrMulPerm;
  
  /* Helper computation variables */
  gsl_matrix *myTmpGradR;
  gsl_matrix *myTmpGradR2;
  gsl_matrix *myTmpR;  
  gsl_vector *myTmpYr;  
  gsl_vector *myTmpCorr;  
  
  const gsl_vector *myP;
  double myPNorm;
  bool isGCD;


  /* Jacobian computation */
  gsl_vector *myTmpJacobianCol;  

  gsl_matrix *myTmpGrad;  

public:

  slraCostFunction( slraStructure *s, int r, const gsl_vector *p, 
                    opt_and_info *opt, gsl_matrix *perm  );
  virtual ~slraCostFunction();
  
  int getD() { return myStruct->getNplusD() - myRank; }
  int getNplusD() { return myStruct->getNplusD(); }
  int getN() { return myRank; }
  int getM() { return myStruct->getM(); }
  int getNp() { return myStruct->getNp(); }

  const gsl_matrix * getPerm() { return myPerm; }
  const gsl_matrix * getSMatr() { return myMatr; }
  
  
  void computeR( gsl_matrix_const_view x_mat, gsl_matrix *R ); 
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
  void computeSr( gsl_matrix *R, gsl_vector *Sr );


  void computeRGammaSr( const gsl_vector *x, gsl_matrix *R, gsl_vector *Sr ) {
    computeR(x, myTmpR);
    myGam->computeCholeskyOfGamma(myTmpR);
    computeSr(myTmpR, Sr);
  } 

  void computeJacobianZij( gsl_vector *res, int i, int j,
                           gsl_vector* yr, gsl_matrix *R, double factor = 0.5 );

  void computePseudoJacobianCorrectFromYr( gsl_vector* yr, gsl_matrix *R, gsl_matrix *jac );
  
  void computePseudoJacobianLsFromYr( gsl_vector* yr, gsl_matrix *R, gsl_matrix *jac );
  void computeGradFromYr( gsl_vector* yr, gsl_matrix *R, gsl_vector *grad );



  void computeCorrectionAndJacobian( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac  );
  void computeFuncAndPseudoJacobianLs( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac );
  void computeFuncAndGrad( const gsl_vector* x, double * f, gsl_vector *grad );
  
  void  computeCorrection( gsl_vector* p, const gsl_vector* x );
  
  static int slra_f_ls( const gsl_vector* x, void* params, gsl_vector *res ) {
    ((slraCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int slra_df_ls( const gsl_vector* x,  void* params, gsl_matrix *jac ) {
    ((slraCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int slra_fdf_ls( const gsl_vector* x,  void* params, gsl_vector *res, gsl_matrix *jac ) {
    ((slraCostFunction *)params)->computeFuncAndPseudoJacobianLs(x, res, jac);
    return GSL_SUCCESS;
  }


  static int slra_f_cor( const gsl_vector* x, void* params, gsl_vector *res ) {
    ((slraCostFunction *)params)->computeCorrectionAndJacobian(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int slra_df_cor( const gsl_vector* x,  void* params, gsl_matrix *jac ) {
    ((slraCostFunction *)params)->computeCorrectionAndJacobian(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int slra_fdf_cor( const gsl_vector* x,  void* params, gsl_vector *res, gsl_matrix *jac ) {
    ((slraCostFunction *)params)->computeCorrectionAndJacobian(x, res, jac);
    return GSL_SUCCESS;
  }


  static double slra_f( const gsl_vector* x, void* params ) {
    double f;
    ((slraCostFunction *)params)->computeFuncAndGrad(x, &f, NULL);
    return f;
  }
  static void  slra_df( const gsl_vector* x, void* params, gsl_vector *grad ) {
    ((slraCostFunction *)params)->computeFuncAndGrad(x, NULL, grad);
  }
  static void  slra_fdf( const gsl_vector* x, void* params, double *f, gsl_vector *grad ) {
    ((slraCostFunction *)params)->computeFuncAndGrad(x, f, grad);
  }
};
