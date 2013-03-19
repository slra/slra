/** Abstract class for s-dependent structure.
 * Assumes that the rows of the matrix are s-dependent, i.e.
 * \f${\bf cov}(S(\widetilde{p})_i, S(\widetilde{p})_j) =const\cdot W_{i,j}\f$
 * and \f$W_{i,j} = 0\f$ for \f$|i-j| < s\f$,
 * where \f$\widetilde{p}\f$ is an uncorrelated sequence with 
 * \f${\bf D} \widetilde{p}_k = \gamma_k\f$.
 *
 * In particular, this implies that \f$\Gamma(R)\f$ is a block matrix with
 * blocks \f$R^{\rm T} W_{i,j} R \f$.
 */
class SDependentStructure : public Structure {
public:
  /** Returns \f$s\f$. */
  virtual size_t getS() const = 0; 
  
  /** Returns \f$res \leftarrow W_{i,j} B\f$ */
  virtual void WijB( gsl_matrix *res, long i, long j, 
                     const gsl_matrix *B ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{i,j} B\f$ */
  virtual void AtWijB( gsl_matrix *res, long i, long j, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWjiB, double beta = 0 ) const = 0;
  /** Returns \f$res \leftarrow \beta res + A^{\rm T} W_{i,j} V\f$ */
  virtual void AtWijV( gsl_vector *res, long i, long j,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWijV, double beta = 0 ) const = 0;
};

