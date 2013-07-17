/** Abstract class for stationary s-dependent structure.
 * A subclass of s-dependent structure, where
 * \f$W_{i,j} = W_{j-i}\f$.
 *
 * In other words, sequence of rows \f$S(\widetilde{p})_j\f$ is
 * stationary.
 */
class StationaryStructure : public SDependentStructure {

public:
  /** Returns \f$res \leftarrow W_{k} B\f$ */
  virtual void WkB( gsl_matrix *res, long k, const gsl_matrix *B ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{k} B\f$ */
  virtual void AtWkB( gsl_matrix *res, long k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWkB, double beta = 0 ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{k} B\f$ */
  virtual void AtWkV( gsl_vector *res, long k,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWkV, double beta = 0 ) const = 0;
                      
  virtual void WijB( gsl_matrix *res, long i, long j, 
                     const gsl_matrix *B ) const {
    WkB(res, j- i, B);
  }
  virtual void AtWijB( gsl_matrix *res, long i, long j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWijB, double beta = 0 ) const {
    AtWkB(res, j - i, A, B, tmpWijB, beta);
  }

  virtual void AtWijV( gsl_vector *res, long i, long j, 
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const {
    AtWkV(res, j - i, A, V, tmpWijV, beta);
  }                      
};
