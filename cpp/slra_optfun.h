class OptFunctionSLRA : public OptFunction {
protected:
  CostFunction &myFun;
  gsl_matrix *myTmpR;  
  gsl_matrix *myPsi;
  gsl_matrix *myTmpXId;
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
public: 
  OptFunctionSLRA( CostFunction &fun, gsl_matrix *Psi );
  virtual ~OptFunctionSLRA();
  virtual int getNvar(); 
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad );

  void computePhat( gsl_vector* p, const gsl_vector* x );

  
  int getRank() { return myPsi->size2 - myFun.getD(); }
  void computeDefaultx( gsl_vector *x ); 

  void RTheta2x( gsl_matrix *RTheta, gsl_vector *x ); 
  void x2RTheta( gsl_matrix *RTheta, const gsl_vector *x ); 

  gsl_matrix x2xmat( const gsl_vector *x );
  static void X2XId( const gsl_matrix *x, gsl_matrix *XId );
  static void PQ2XId( const gsl_matrix *PQ, gsl_matrix * x );
};


class OptFunctionSLRACholesky : public OptFunctionSLRA {
public: 
  OptFunctionSLRACholesky( CostFunction &fun, gsl_matrix *psi ) : 
      OptFunctionSLRA(fun, psi) {}
  virtual ~OptFunctionSLRACholesky() {}
  virtual int getNsq() { return myFun.getN() * myFun.getD(); }
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac );
};

class OptFunctionSLRACorrection : public OptFunctionSLRA {
public: 
  OptFunctionSLRACorrection( CostFunction &fun, gsl_matrix *psi ) : 
      OptFunctionSLRA(fun, psi)  {}
  virtual ~OptFunctionSLRACorrection() {}
  virtual int getNsq() { return myFun.getNp(); }
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac );
};

