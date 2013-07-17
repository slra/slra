/** Implementation of Cholesky for StripedStructure. */
class StripedCholesky : virtual public Cholesky {
  Cholesky **myGamma;
  size_t myD;
  size_t myNGamma;
  const StripedStructure *myS;
public:  
  /** Constructs a stripe of Cholesky objects
   * using createDGamma  for each block of the stripe */
  StripedCholesky( const StripedStructure *s, size_t D );
  virtual ~StripedCholesky();

  /** @name Implementing Cholesky interface */
  /**@{*/  
  virtual void calcGammaCholesky( const gsl_matrix *R, double reg_gamma );
  virtual void multInvCholeskyVector( gsl_vector * yr, long trans );  
  virtual void multInvGammaVector( gsl_vector * yr );                
  /**@}*/
};

