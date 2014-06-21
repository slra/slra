/** Abstract class for stationary \f$\mu\f$-dependent Structure.
 * It a pair (structure, weights) such that the matrix \f$\mathrm{V}\f$
 * (see eqn. \f$(\mathrm{V})\f$ in \cite slra-efficient) 
 * is block-Toeplitz block-banded, i.e.
 * \f$\mathrm{V}_{\#ij} = \mathrm{V}_{j-i}\f$. The matrices
 * \f$\mathrm{V}_{k}\f$ are defined for \f$-\mu \le k \le \mu\f$.
 *
 * A special case of StationaryStructure is considered 
 * in (b) of Theorem 1 in \cite slra-efficient. 
 */
class StationaryStructure : public MuDependentStructure {

public:
  /** Returns \f$X \leftarrow \mathrm{V}_{k} B\f$ 
   * for \f$B \in \mathbb{R}^{m \times Q}\f$. */
  virtual void VkB( gsl_matrix *X, ///< [out] the \f$m \times Q\f$ matrix \f$X\f$
                    long k,        ///< [in] index of block diagonal \f$k\f$ 
                    const gsl_matrix *B ///< [in] matrix \f$B\f$
                   ) const = 0;

  /** Updates \f$X \leftarrow \beta X + A^{\top} \mathrm{V}_{k} B\f$, 
   * for \f$A \in \mathbb{R}^{m \times P}\f$, \f$B \in \mathbb{R}^{m \times Q}\f$. */ 
  virtual void AtVkB( gsl_matrix *X,  ///< [out,in] the \f$P \times Q\f$ matrix \f$X\f$ 
                      long k,        ///< [in] index of block diagonal \f$k\f$ 
                      const gsl_matrix *A, ///< [in]  matrix \f$A\f$
                      const gsl_matrix *B, ///< [in]  matrix \f$B\f$
                      gsl_matrix *tmpVkB,  ///< [out] temporary \f$m \times Q\f$  matrix 
                      double beta = 0      ///< scalar \f$\beta\f$
                     ) const = 0;

  /** Updates \f$u \leftarrow \beta u + A^{\top} \mathrm{V}_{k} v\f$, 
   * for \f$A \in \mathbb{R}^{m \times P}\f$. */
  virtual void AtVkV( gsl_vector *u, ///< [out,in] vector \f$u \in \mathbb{R}^P\f$ 
                      long k,        ///< [in] index of block diagonal \f$k\f$ 
                      const gsl_matrix *A, ///< [in]  matrix \f$A\f$
                      const gsl_vector *v, ///< [in] vector \f$v \in \mathbb{R}^m\f$ 
                      gsl_vector *tmpVkV,  ///< [out] temporary vector of length \f$m\f$ 
                      double beta = 0      ///< scalar \f$\beta\f$
                     ) const = 0;

  virtual void VijB( gsl_matrix *X, long i_1, long j_1, 
                     const gsl_matrix *B ) const {
    VkB(X, j_1 - i_1, B);
  }
  virtual void AtVijB( gsl_matrix *X, long i_1, long j_1, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpVijB, double beta = 0 ) const {
    AtVkB(X, j_1 - i_1, A, B, tmpVijB, beta);
  }

  virtual void AtVijV( gsl_vector *u, long i_1, long j_1, 
                      const gsl_matrix *A, const gsl_vector *v, 
                      gsl_vector *tmpVijV, double beta = 0 ) const {
    AtVkV(u, j_1 - i_1, A, v, tmpVijV, beta);
  }                      
};
