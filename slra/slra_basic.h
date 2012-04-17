class slraException {
  static const int MSG_MAX = 200;

  char myMsg[MSG_MAX];
public:
  slraException( const char *msg, ... );
  const char *getMessage()  { return myMsg; }
};


class slraGammaCholesky {
public:  
  virtual  ~slraGammaCholesky() {}
  virtual void computeCholeskyOfGamma( gsl_matrix *R ) = 0;

  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans ) = 0;  
  virtual void multiplyInvGammaVector( gsl_vector * yr ) = 0;                
  virtual void multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans ) { /* Default implementation */
    for (int i = 0; i < yr_matr->size1; i++) {
      gsl_vector_view row = gsl_matrix_row(yr_matr, i);
      multiplyInvCholeskyVector(&row.vector, trans);
    }
  }
};

class slraDGamma {
public:  
  virtual ~slraDGamma() {}
  virtual void computeYrtDgammaYr( gsl_matrix *grad, gsl_matrix *R, gsl_vector *yr ) = 0;

  virtual void computeDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   gsl_matrix *perm, int i, int j, gsl_vector *Yr ) = 0;

};


class slraStructure {
public:
  virtual ~slraStructure() {}
  virtual int getNp() const = 0;
  virtual int getNplusD() const = 0;
  virtual int getM() const = 0;
  virtual int getS() const = 0;
  
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )  = 0; 
  
  virtual slraGammaCholesky *createGammaComputations( int r, double reg_gamma ) const = 0;
  virtual slraDGamma *createDerivativeComputations( int r ) const = 0;
  virtual void correctVector( gsl_vector* p, gsl_matrix *R, gsl_vector *yr ) = 0;
};

class slraStationaryStructure : public slraStructure {
public:
  virtual const gsl_matrix *getWk( int k )  const = 0; 
};


