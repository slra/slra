/** Implementation of DGamma for StripedStructure */
class StripedDGamma : virtual public DGamma {
  DGamma **myLHDGamma;
  const StripedStructure *myS;
  gsl_matrix *myTmpGrad;
public:  
  /** Constructs a stripe of DGamma objects
   * using createDGamma  for each block of the stripe */
  StripedDGamma( const StripedStructure *s, size_t D  ) ;
  virtual ~StripedDGamma();

  /** @name Implementing DGamma interface */
  /**@{*/
  virtual void calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                   const gsl_vector *yr );
  virtual void calcDijGammaYr( gsl_vector *res, const gsl_matrix *R, 
                   size_t i, size_t j, const gsl_vector *Yr );
  /**@}*/
};

