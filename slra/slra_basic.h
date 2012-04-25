class Exception {
  static const int MSG_MAX = 200;

  char myMsg[MSG_MAX];
public:
  Exception( const char *msg, ... );
  const char *getMessage()  { return myMsg; }
};


class Cholesky {
public:  
  virtual  ~Cholesky() {}
  virtual void calcGammaCholesky( gsl_matrix *R ) = 0;

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans ) = 0;  
  virtual void multiplyInvGammaVector( gsl_vector * yr ) = 0;                
  virtual void multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans );
};

class DGamma {
public:  
  virtual ~DGamma() {}
  virtual void calcYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, 
                   gsl_vector *yr ) = 0;
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr ) = 0;
};


class Structure {
public:
  virtual ~Structure() {}
  virtual int getNp() const = 0;
  virtual int getNplusD() const = 0;
  virtual int getM() const = 0;
  
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  = 0; 
  
  virtual Cholesky *createCholesky( int D, double reg_gamma ) const = 0;
  virtual DGamma *createDGamma( int D ) const = 0;
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) = 0;
};



class SDependentStructure : public Structure {
public:
  virtual int getS() const = 0;
  virtual void WijB( gsl_matrix *res, int i, int j, const gsl_matrix *B ) const = 0;
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const = 0;
  virtual void AtWijV( gsl_vector *res, int i, int j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const = 0;
};


class StationaryStructure : public SDependentStructure {

public:
  virtual void WkB( gsl_matrix *res, int k, const gsl_matrix *B ) const = 0;
  virtual void AtWkB( gsl_matrix *res, int k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWkB, double beta = 0 ) const = 0;
  virtual void AtWkV( gsl_vector *res, int k,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWkV, double beta = 0 ) const = 0;
                      
  virtual void WijB( gsl_matrix *res, int i, int j, const gsl_matrix *B ) const {
    WkB(res, j- i, B);
  }
  virtual void AtWijB( gsl_matrix *res, int i, int j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWijB, double beta = 0 ) const {
    AtWkB(res, j - i, A, B, tmpWijB, beta);
  }

  virtual void AtWijV( gsl_vector *res, int i, int j, 
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const {
    AtWkV(res, j - i, A, V, tmpWijV, beta);
  }                      
};
