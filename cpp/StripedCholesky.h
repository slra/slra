/** Implementation of Cholesky for StripedStructure. 
 */
class StripedCholesky : virtual public Cholesky {
  Cholesky **myGamma;
  size_t myD;
  size_t myNGamma;
  const StripedStructure *myStruct;
public:  
  /** Constructs a StripedCholesky from StripedStructure.
   * Create a stripe of Cholesky objects
   * using createCholesky  for each block of the stripe. 
   * Saves memory if StripedStructure::isSameGamma() returns `true`.
   * @param[in] s    Pointer to the corresponding StripedStructure.
   * @param[in] d    number of rows \f$d\f$ of  the matrix \f$R\f$ */
  StripedCholesky( const StripedStructure *s, size_t d );
  virtual ~StripedCholesky();

  /** @name Implementing Cholesky interface */
  /**@{*/  
  virtual void calcGammaCholesky( const gsl_matrix *Rt, double reg );
  virtual void multInvCholeskyVector( gsl_vector * yr, long trans );  
  virtual void multInvGammaVector( gsl_vector * yr );                
  /**@}*/
};

