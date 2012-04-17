class slraGammaCholeskyBTBanded : public slraGammaCholesky {
protected:
  int my_use_slicot;
  double my_reg_gamma;
  
  const slraWkInterface *myW;
  
  size_t myN, myD;
  
  size_t myMg;
  
  size_t s_minus_1;
  size_t d_times_s;
  size_t d_times_Mg;
  size_t d_times_s_minus_1;
  
  int myCholeskyWorkSize;
 
  double *myGammaVec;
  gsl_matrix *myGamma;
  gsl_matrix *myWkTmp;
  double *myPackedCholesky;
  double *myCholeskyWork;
  
public:
  slraGammaCholeskyBTBanded( const slraWkInterface *s, int r, int Mg,
     int use_slicot, double reg_gamma  );
  virtual ~slraGammaCholeskyBTBanded();

  int getD() const { return myD; }

  virtual void computeCholeskyOfGamma( gsl_matrix *R );

  virtual void multiplyInvPartCholeskyArray( double * yr, int trans, size_t size, size_t chol_size );
  virtual void multiplyInvPartGammaArray( double * yr, size_t size, size_t chol_size );


  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans ) {
    if (yr->stride != 1) {
      throw new slraException("Cannot multiply vectors with stride != 1\n");
    }
    multiplyInvPartCholeskyArray(yr->data, trans, yr->size, d_times_Mg);
  }
  virtual void multiplyInvGammaVector( gsl_vector * yr ) {
    if (yr->stride != 1) {
      throw new slraException("Cannot multiply vectors with stride != 1\n");
    }
    multiplyInvPartGammaArray(yr->data, yr->size, d_times_Mg);
  }

  virtual void multiplyInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans ) {
    if (yr_matr->size2 != yr_matr->tda) {
      slraGammaCholesky::multiplyInvCholeskyTransMatrix(yr_matr, trans);
    } else {
      multiplyInvPartCholeskyArray(yr_matr->data, trans, yr_matr->size1 * yr_matr->size2, d_times_Mg);
    }
  }

};

class slraGammaCholeskySameDiagBTBanded : public slraGammaCholesky {
  slraGammaCholeskyBTBanded *myBase;
  const slraMosaicHankelStructure *myStruct;
public:  
  slraGammaCholeskySameDiagBTBanded( const slraMosaicHankelStructure *s, int r, int use_slicot, double reg_gamma  ) :
      myStruct(s) {
    myBase = (slraGammaCholeskyBTBanded *)myStruct->getMaxBlock()->createGammaComputations(r, reg_gamma);  
  }
  virtual ~slraGammaCholeskySameDiagBTBanded() {
    delete myBase;
  }
  
  virtual void computeCholeskyOfGamma( gsl_matrix *R ) {
    myBase->computeCholeskyOfGamma(R);
  }
  
  virtual void multiplyInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multiplyInvGammaVector( gsl_vector * yr );                
};

