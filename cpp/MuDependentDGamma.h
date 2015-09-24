/** Implementation of DGamma class for the MuDependentStructure. 
 * MuDependentDGamma::calcYrtDgammaYr() is implemented using
 * eqn. \f$(\nabla_{d\times m}\f) in \cite slra-efficient.
 */
class MuDependentDGamma : public DGamma {
private:
  const MuDependentStructure *myStruct;
  size_t myD;
  gsl_vector *myTmp1, *myTmp2, *myTmp3;
  gsl_vector *myYrR;
  gsl_matrix *myEye;
public:  
  /** Constructs a MuDependentDGamma object.
   * @copydetails MuDependentCholesky::MuDependentCholesky */
  MuDependentDGamma( const MuDependentStructure *s, size_t d );
  virtual ~MuDependentDGamma();

  /** @name Implementing DGamma interface */
  /**@{*/  
  virtual void calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt, 
                              const gsl_matrix *Yt );
  virtual void calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt, 
                   size_t j_1, size_t i_1, const gsl_vector *y,
                   const gsl_matrix *Phi = NULL );
  /**@}*/  
};
