class OptFunctionSLRA : public OptFunction {
protected:
  CostFunction &myFun;
  gsl_matrix *myTmpR;  
  gsl_matrix *myPerm;
  gsl_matrix *myPsi, *myPhi;
  gsl_matrix *myTmpXId;
  void computeR( const gsl_vector *x, gsl_matrix *R ); 
public: 
  OptFunctionSLRA( CostFunction &fun, gsl_matrix *Phi, gsl_matrix *Psi );
  virtual ~OptFunctionSLRA();
  virtual int getNvar(); 
  virtual void computeFuncAndGrad( const gsl_vector* x, double* f, gsl_vector *grad );

  void computeCorrection( gsl_vector* p, const gsl_vector* x );

  
  int getRank() { return myPerm->size2 - myFun.getD(); }
  gsl_matrix * getPerm() { return myPerm; }
  void computeDefaultRTheta( gsl_matrix *RTheta ); 
  void computeDefaultx( gsl_vector *x ); 

  void RTheta2x( gsl_matrix *RTheta, gsl_vector *x ); 
  void x2RTheta( gsl_matrix *RTheta, const gsl_vector *x ); 

  static void X2XId( const gsl_matrix *x, gsl_matrix *XId );
  static void PQ2XId( const gsl_matrix *PQ, gsl_matrix * x );
};


class OptFunctionSLRACholesky : public OptFunctionSLRA {
public: 
  OptFunctionSLRACholesky( CostFunction &fun, gsl_matrix *phi, gsl_matrix *psi ) : 
      OptFunctionSLRA(fun, phi, psi) {}
  virtual ~OptFunctionSLRACholesky() {}
  virtual int getNsq() { return myFun.getSt()->getN() * myFun.getD(); }
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac );
};

class OptFunctionSLRACorrection : public OptFunctionSLRA {
public: 
  OptFunctionSLRACorrection( CostFunction &fun, gsl_matrix *phi, gsl_matrix *psi ) : 
      OptFunctionSLRA(fun, phi, psi)  {}
  virtual ~OptFunctionSLRACorrection() {}
  virtual int getNsq() { return myFun.getSt()->getNp(); }
  virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res, gsl_matrix *jac );
};

