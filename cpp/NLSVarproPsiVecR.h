class NLSVarproPsiVecR : public NLSVarpro {
protected:
  size_t myNEssVar;
  gsl_matrix *myPsiTBig;
  gsl_vector *myTmpRVec;
  gsl_matrix myTmpR;
public:
  NLSVarproPsiVecR( VarproFunction &fun,
                    const gsl_matrix *PsiT );
  virtual ~NLSVarproPsiVecR();
  
  virtual size_t getNvar() { return myPsiTBig->size2; }
  virtual size_t getNEssVar() { return myNEssVar; }

  virtual void RTheta2x( const gsl_matrix *RTheta, gsl_vector *x );
  
  virtual void x2RTheta( gsl_matrix *RTheta, const gsl_vector *x );
  
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad ) {
    double tmp;
    
    
    x2RTheta(&myTmpR, x);
    if (grad == NULL) {
      myFun.computeFuncAndGrad(&myTmpR, f, NULL, NULL);
    } else {
      gsl_matrix gradV = gsl_matrix_view_vector(grad, 1, grad->size).matrix;
      myFun.computeFuncAndGrad(&myTmpR, f, myPsiTBig, &gradV);
    }
  }
};

class NLSVarproPsiVecRCholesky : public NLSVarproPsiVecR {
public:
    NLSVarproPsiVecRCholesky( VarproFunction &fun, const gsl_matrix *PsiT ) :  NLSVarproPsiVecR(fun, PsiT) {}
    virtual ~NLSVarproPsiVecRCholesky() {}
    virtual size_t getNsq() { return myFun.getN() * myFun.getD(); }
    virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res,
                                   gsl_matrix *jac ) {
      x2RTheta(&myTmpR, x);
      myFun.computeFuncAndPseudoJacobianLs(&myTmpR, myPsiTBig, res, jac);
    }
};

class NLSVarproPsiVecRCorrection : public NLSVarproPsiVecR {
public:
    NLSVarproPsiVecRCorrection( VarproFunction &fun, const gsl_matrix *PsiT ) : NLSVarproPsiVecR(fun, PsiT)  {}
    virtual ~NLSVarproPsiVecRCorrection() {}
    virtual size_t getNsq() { return myFun.getNp(); }
    virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res,
                                   gsl_matrix *jac ) {
      x2RTheta(&myTmpR, x);
      myFun.computeCorrectionAndJacobian(&myTmpR, myPsiTBig, res, jac);
    }
};


