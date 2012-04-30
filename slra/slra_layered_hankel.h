class LayeredHStructure : public StationaryStructure {
  typedef struct {
    size_t blocks_in_row;       /* Number of blocks in a row of Ci */
    double inv_w;            /* Square root of inverse of the weight */
  } Layer;

  int myQ;	                /* number of layers */
  size_t myM;
  size_t myNplusD;
  size_t myMaxLag;
  gsl_matrix **myA;
  Layer *mySA;	/* q-element array describing C1,...,Cq; */  

  void computeStats();
  void computeWkParams(); 
protected:
  int nvGetNp() const { return (myM - 1) * myQ + myNplusD; }  
public:
  LayeredHStructure( const double *oldNk, size_t q, int M, 
                     const double *layer_w = NULL );
  virtual ~LayeredHStructure();

  /* Structure methods */
  virtual int getNplusD() const { return myNplusD; }
  virtual int getM() const { return myM; }
  virtual int getNp() const { return nvGetNp(); }
  virtual Cholesky *createCholesky( int D, double reg_gamma ) const;
  virtual DGamma *createDGamma( int D ) const;
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ); 
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr );

  /* StationaryStructure methods */
  virtual int getS() const { return myMaxLag; }
  virtual void WkB( gsl_matrix *res, int k, const gsl_matrix *B ) const;
  virtual void AtWkB( gsl_matrix *res, int k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWkB, double beta = 0 ) const;
  virtual void AtWkV( gsl_vector *res, int k,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWkV, double beta = 0 ) const;

  /* Structure-specific methods */
  const gsl_matrix *getWk( int l ) const { return myA[l]; }
  void setM( int m );
  int getQ() const { return myQ; }
  int getMaxLag() const { return myMaxLag; }
  int getLayerLag( int l ) const { return mySA[l].blocks_in_row; }
  bool isLayerExact( int l ) const { return (mySA[l].inv_w == 0.0); }
  double getLayerInvWeight( int l ) const { return mySA[l].inv_w; }
  int getLayerNp( int l ) const { return getLayerLag(l) + getM() - 1; }
};

class MosaicHStructure : public StripedStructure {
  bool myWkIsCol;
protected:
  static Structure **allocStripe( gsl_vector *oldNk, gsl_vector *oldMl,  
                                  gsl_vector *Wk, bool wkIsCol = false );
public:
  MosaicHStructure( gsl_vector *oldNk, gsl_vector *oldMl,  
                    gsl_vector *Wk, bool wkIsCol = false );
  virtual ~MosaicHStructure() {}
  virtual Cholesky *createCholesky( int r, double reg_gamma ) const;
};


