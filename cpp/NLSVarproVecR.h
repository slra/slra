class NLSVarproVecR : public NLSVarpro {
protected:
  VarproFunction &myFun;
public: 
  NLSVarproVecR( VarproFunction &fun ) : myFun(fun) {}
  virtual ~NLSVarproVecR() {}
  
  virtual size_t getD() { return myFun.getD(); }
  virtual size_t getM() { return myFun.getNrow(); }

  virtual size_t getNsq() { return myFun.getN() * myFun.getD(); }
  virtual size_t getNvar() { return myFun.getNrow() * myFun.getD(); }
  virtual size_t getNEssVar() { return (myFun.getNrow() - myFun.getD()) * myFun.getD(); }

  virtual gsl_matrix x2xmat( const gsl_vector *x ) {
    return gsl_matrix_const_view_vector(x, getM(), getD()).matrix;
  }

  virtual void computePhat( gsl_vector* p, const gsl_vector* x ) {
    gsl_matrix x_mat = x2xmat(x);
	myFun.computePhat(p, &x_mat);
  }
  
  virtual void computeDefaultx( gsl_vector *x ) {
    gsl_matrix x_mat = x2xmat(x);
    myFun.computeDefaultRTheta(&x_mat);
  }


  
  virtual void RTheta2x( gsl_matrix *RTheta, gsl_vector *x ) {
    gsl_matrix x_mat = x2xmat(x);
    gsl_matrix_memcpy(&x_mat, RTheta);
  }
  
  virtual void x2RTheta( gsl_matrix *RTheta, const gsl_vector *x ) {
    gsl_matrix x_mat = x2xmat(x);
    gsl_matrix_memcpy(RTheta, &x_mat);
  }
  
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad ) {
    gsl_matrix tmpR = x2xmat(x);										 
    if (grad == NULL) {
      myFun.computeFuncAndGrad(&tmpR, f, NULL, NULL);
    } else {
      gsl_matrix gradM = x2xmat(grad);
      myFun.computeFuncAndGrad(&tmpR, f, NULL, &gradM);
    }
  }


  
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, 
                                   gsl_matrix *jac ) {
    gsl_matrix tmpR = x2xmat(x);
    myFun.computeFuncAndPseudoJacobianLs(&tmpR, NULL, res, jac); 
  }   
};
