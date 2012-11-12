/** Prototype class for striped structure.
 * Abstract class that facilitates creation of structures of form
 * \f$\mathcal{S}(p) = 
 * \begin{bmatrix} 
 *   \mathcal{S}_1 (p^{(1)})\\ \vdots \\\mathcal{S}_N (p^{(N)})
 * \end{bmatrix}\f$
 * All structures must have the same \f$m\f$.
 */
class StripedStructure : public Structure {
  size_t myBlocksN;
  Structure **myStripe;

  /* Helper variables */
  size_t myN;
  size_t myNp;
  size_t myMaxNkInd;
protected:
  /** Constructs striped structure from array.
   * @param blocksN \f$N\f$ 
   * @param stripe array of Structure objects
   */
  StripedStructure( size_t blocksN, Structure *stripe[] );
public:
  virtual ~StripedStructure();

  /** @name Implementing Structure interface */
  /**@{*/
  virtual int getN() const { return myN; }
  virtual int getNp() const { return myNp; }
  virtual int getM() const { return myStripe[0]->getM(); }
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) ;
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr,
                         int wdeg = 2 );
  virtual Cholesky *createCholesky( int D, double reg_gamma ) const;
  virtual DGamma *createDGamma( int D ) const;
  /**@}*/
  
  /** @name Structure-specific methods */
  /**@{*/
  int getBlocksN() const { return myBlocksN; } /**< Returns \f$N\f$ */
  /** Returns \f$\mathcal{S}_{k}\f$ */
  const Structure *getBlock( size_t k ) const { 
    return myStripe[k]; 
  }
  /** Returns \f$\mathcal{S}_{k_{max}}\f$, where 
   * \f$k_{max} = argmax \{n_k\}_{k=1}^{N}\f$  */
  const Structure *getMaxBlock() const { 
    return getBlock(myMaxNkInd); 
  }
  /**@}*/
};


/** Implementation of Cholesky for StripedStructure. */
class StripedCholesky : virtual public Cholesky {
  Cholesky **myGamma;
  int myD;
  const StripedStructure *myS;
public:  
  /** Constructs a stripe of Cholesky objects
   * using createDGamma  for each block of the stripe */
  StripedCholesky( const StripedStructure *s, int D, double reg_gamma );
  virtual ~StripedCholesky();

  /** @name Implementing Cholesky interface */
  /**@{*/  
  virtual void calcGammaCholesky( const gsl_matrix *R, bool regularize = true );
  virtual void multInvCholeskyVector( gsl_vector * yr, int trans );  
  virtual void multInvGammaVector( gsl_vector * yr );                
  /**@}*/
};

/** Implementation of DGamma for StripedStructure */
class StripedDGamma : virtual public DGamma {
  DGamma **myLHDGamma;
  const StripedStructure *myS;
  gsl_matrix *myTmpGrad;
public:  
  /** Constructs a stripe of DGamma objects
   * using createDGamma  for each block of the stripe */
  StripedDGamma( const StripedStructure *s, int D  ) ;
  virtual ~StripedDGamma();

  /** @name Implementing DGamma interface */
  /**@{*/
  virtual void calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                   const gsl_vector *yr );
  virtual void calcDijGammaYr( gsl_vector *res, gsl_matrix *R, 
                   int i, int j, gsl_vector *Yr );
  /**@}*/
};

