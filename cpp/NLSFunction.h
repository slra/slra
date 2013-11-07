class OptFunction {
public:
  virtual ~OptFunction() {}
  virtual size_t getNvar() = 0;
  virtual size_t getNsq() = 0;
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, 
                                   gsl_vector *grad ) = 0;
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, 
                                  gsl_matrix *jac ) = 0;

  static double _f( const gsl_vector* x, void* params ) {
    double f;
    ((OptFunction *)params)->computeFuncAndGrad(x, &f, NULL);
    return f;
  }
  static void _df( const gsl_vector* x, void* params, gsl_vector *grad ) {
    ((OptFunction *)params)->computeFuncAndGrad(x, NULL, grad);
  }
  static void _fdf( const gsl_vector* x, void* params, double* f, 
                     gsl_vector* grad ) {
    ((OptFunction *)params)->computeFuncAndGrad(x, f, grad);
  }

  static int _f_ls( const gsl_vector* x, void* params, gsl_vector* res ) {
    ((OptFunction *)params)->computeFuncAndJac(x, res, NULL);
    return GSL_SUCCESS;
  }
  static int _df_ls( const gsl_vector* x,  void* params, gsl_matrix* jac ) {
    ((OptFunction *)params)->computeFuncAndJac(x, NULL, jac);
    return GSL_SUCCESS;
  }
  static int _fdf_ls( const gsl_vector* x, void* params, gsl_vector*res, 
                       gsl_matrix *jac ) {
    ((OptFunction *)params)->computeFuncAndJac(x, res, jac);
    return GSL_SUCCESS;
  }
};



