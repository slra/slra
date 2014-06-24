/** Implementation of DGamma for StripedStructure */
class StripedDGamma : virtual public DGamma {
  DGamma **myLHDGamma;
  const StripedStructure *myS;
  gsl_matrix *myTmpGrad;
public:  
  /** Constructs a stripe of DGamma objects
   * using createDGamma  for each block of the stripe. */
  StripedDGamma( const StripedStructure *s, size_t d  ) ;
  virtual ~StripedDGamma();

  /** @name Implementing DGamma interface */
  /**@{*/
  virtual void calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt, 
                   const gsl_matrix *Yt );
  virtual void calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt, 
                    size_t j_1, size_t i_1, const gsl_vector *y );
  /**@}*/
};

