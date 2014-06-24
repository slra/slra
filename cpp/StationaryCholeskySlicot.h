#ifdef USE_SLICOT
/** Implementation of Cholesky class for the StationaryStructure
 * using the function MB02GD of the SLICOT library. 
 *
 * This class is available only if the macro USE_SLICOT is defined. 
 */
class StationaryCholeskySlicot : public StationaryCholesky {
public:
  /** Constructs the StationaryCholeskySlicot object.
   * @copydetails  StationaryCholesky::StationaryCholesky */
  StationaryCholeskySlicot( const StationaryStructure *s, size_t d );
  virtual ~StationaryCholeskySlicot();
  
  /** @name Implementing Cholesky interface */
  /**@{*/
  virtual void calcGammaCholesky( const gsl_matrix *Rt, double reg = 0 );
  /**@}*/
private:
  double *myGammaVec;
  double *myCholeskyWork;
  size_t myCholeskyWorkSize;
};
#endif /* USE_SLICOT */

