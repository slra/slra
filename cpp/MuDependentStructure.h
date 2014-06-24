/** Abstract class for \f$\mu\f$-dependent Structure.
 * The Structure (structure and weight pair) is called \f$\mu\f$-dependent if
 * and \f$\mathrm{V}_{\#ij} = 0\f$ (or equivalently \f$\mathrm{\Gamma}_{\#ij} = 0\f$),
 * for all \f$i,j: |i-j| \ge \mu\f$. (See also eqs. \f$(\mathrm{V})\f$ and
 * \f$(\Gamma_{\#ij})\f$ in \cite slra-efficient.)
 *
 * Thus \f$\mu\f$ is the block bandwidth of \f$\Gamma(R)\f$.
 * See also Theorems 1--3 in \cite slra-efficient.
 *
 * The notion of \f$\mu\f$-dependence also coincides with the notion of
 * \f$s\f$-dependence from \cite KMV02e.
 */
class MuDependentStructure : public Structure {
public:
  /** Returns \f$\mu\f$. */
  virtual size_t getMu() const = 0; 
  
  /** Returns \f$X \leftarrow \mathrm{V}_{\#ij} B\f$ 
   * for \f$B \in \mathbb{R}^{m \times Q}\f$. */
  virtual void VijB( gsl_matrix *X, ///< [out] the \f$m \times Q\f$ matrix \f$X\f$
                     long i_1,   ///< [in] \f$0\f$-based index \f$i_1\f$, such that 
                                 ///   \f$i=i_1+1\f$  and \f$0 \le i_1 <n \f$.   
                     long j_1,   ///< [in] \f$0\f$-based index \f$j_1\f$, such that 
                                 ///   \f$j=j_1+1\f$  and \f$0 \le j_1 <n \f$.   
                     const gsl_matrix *B ///< [in] matrix \f$B\f$
                    ) const = 0;

  /** Updates \f$X \leftarrow \beta X + A^{\top} \mathrm{V}_{\#ij} B\f$, 
   * for \f$A \in \mathbb{R}^{m \times P}\f$, \f$B \in \mathbb{R}^{m \times Q}\f$. */ 
  virtual void AtVijB( gsl_matrix *X,  ///< [out,in] the \f$P \times Q\f$ matrix \f$X\f$ 
                       long i_1,   ///< [in] \f$0\f$-based index \f$i_1\f$, such that 
                                   ///   \f$i=i_1+1\f$  and \f$0 \le i_1 <n \f$.   
                       long j_1,   ///< [in] \f$0\f$-based index \f$j_1\f$, such that 
                                   ///   \f$j=j_1+1\f$  and \f$0 \le j_1 <n \f$.   
                       const gsl_matrix *A, ///< [in]  matrix \f$A\f$
                       const gsl_matrix *B, ///< [in]  matrix \f$B\f$
                       gsl_matrix *tmpVijB, ///< [out] temporary \f$m \times Q\f$  matrix 
                       double beta = 0      ///< scalar \f$\beta\f$
                      ) const = 0;

  /** Updates \f$u \leftarrow \beta u + A^{\top} \mathrm{V}_{\#ij} v\f$, 
   * for \f$A \in \mathbb{R}^{m \times P}\f$. */
  virtual void AtVijV( gsl_vector *u, ///< [out,in] vector \f$u \in \mathbb{R}^P\f$ 
                       long i_1,   ///< [in] \f$0\f$-based index \f$i_1\f$, such that 
                                   ///   \f$i=i_1+1\f$  and \f$0 \le i_1 <n \f$.   
                       long j_1,   ///< [in] \f$0\f$-based index \f$j_1\f$, such that 
                                   ///   \f$j=j_1+1\f$  and \f$0 \le j_1 <n \f$.   
                       const gsl_matrix *A, ///< [in]  matrix \f$A\f$
                       const gsl_vector *v, ///< [in] vector \f$v \in \mathbb{R}^m\f$ 
                       gsl_vector *tmpVijV, ///< [out] temporary vector of length \f$m\f$ 
                       double beta = 0      ///< scalar \f$\beta\f$
                      ) const = 0;
};

