/** Implementation of DGamma class for the StationaryStructure. 
 * StationaryDGamma::calcYtDgammaY() is implemented using
 * eqn. Proposition 5 in \cite slra-efficient.
 */
class StationaryDGamma : public DGamma {
private:
  const StationaryStructure *myW;
  size_t  myD;
  
  gsl_vector *myTempVkColRow;
  gsl_vector *myDGammaVec;
  gsl_matrix *myDGammaTrMat;
  gsl_matrix *myDGamma;
  gsl_vector *myTmpCol;
  
  gsl_matrix *myVk_R;
  gsl_matrix *myN_k;
  gsl_matrix *myEye;
public:
  StationaryDGamma( const StationaryStructure *s, size_t D );
  virtual ~StationaryDGamma();
  size_t getD() const { return myD; }

  
  /** @name Implementing DGamma interface */
  /**@{*/  
  virtual void calcYtDgammaY(  gsl_matrix *At, const gsl_matrix *Rt, 
                              const gsl_matrix *Yt );
  virtual void calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt, 
                    size_t j_1, size_t i_1, const gsl_vector *y );
  /**@}*/
};



