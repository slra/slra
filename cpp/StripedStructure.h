/** Prototype class for striped structure and block-diagonal weight matrix.
 * Abstract class for creation of Structure object for a striped structure
 * \f[
 *    \mathscr{S}(p) = 
 * \begin{bmatrix} 
 *   \mathscr{S}^{(1)} (p^{(1)}) & \cdots & \mathscr{S}^{(N)} (p^{(N)})
 * \end{bmatrix}
 * \f]
 * and block-diagonal weight matrix
 * \f[\mathrm{W}=\mathrm{blkdiag}(\mathrm{W}^{(1)},\ldots,\mathrm{W}^{(N)}),\f]
 * where \f$\mathscr{S}^{(l)}: \mathbb{R}^{n_p^{(l)}} \to \mathbb{R}^{m\times n_l}\f$ and
 * \f$\mathrm{W}^{(l)} \in \mathbb{R}^{n_p^{(l)} \times n_p^{(l)}}\f$
 * 
 * For more details, see Lemma 1 in \cite slra-efficient.
 */
class StripedStructure : public Structure {
  size_t myBlocksN;
  Structure **myStripe;

  /* Helper variables */
  size_t myN;
  size_t myNp;
  size_t myMaxNlInd;
  bool myIsSameGamma;
public:
  /** Constructs StripedStructure from array of Structure  objects.
   * @param blocksN     \f$N\f$ --- number of blocks
   * @param stripe      array of Structure objects
   * @param isSameGamma if `isSameGamma == true`, it is assumed that all 
   *     \f$\Gamma_{\mathscr{S}^{(l)}}(R)\f$  are submatrices of 
   *     \f$\Gamma_{\mathscr{S}^{(l_{max})}}(R)\f$ (where \f$l_{max}\f$ is defined  
   *     in StripedStructure::getMaxBlock). In this case, memory is saved in 
   *     StripedCholesky and StripedDGamma objects.
   */
  StripedStructure( size_t blocksN, Structure *stripe[], 
                    bool isSameGamma = false );
  virtual ~StripedStructure();

  /** @name Implementing Structure interface */
  /**@{*/
  virtual size_t getN() const { return myN; }
  virtual size_t getNp() const { return myNp; }
  virtual size_t getM() const { return myStripe[0]->getM(); }
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) ;
  
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
                                   const gsl_vector *y, 
                                   double alpha = -1, double beta = 1,
                                   bool skipFixedBlocks = true ); 
  virtual void multByWInv( gsl_vector* p, long deg = 2 );
  virtual Cholesky *createCholesky( size_t d ) const;
  virtual DGamma *createDGamma( size_t d ) const;
  /**@}*/
  
  /** @name StripedStructure-specific methods */
  /**@{*/
  size_t getBlocksN() const { return myBlocksN; } /**< Returns \f$N\f$ */
  
  /** Returns a Structure object for the pair \f$\mathscr{S}^{(l)},\mathrm{W}^{(l)}\f$.
   * @param l_1 index \f$0\f$-based index \f$l_1\f$ such that \f$l = l_1+1\f$ and
   * \f$0 \le l_1 <N \f$.    
   */
  const Structure *getBlock( size_t l_1 ) const { 
    return myStripe[l_1]; 
  }

  /** Returns \f$\mathscr{S}^{(l_{max})}\f$, where 
   * \f$l_{max} := \mathrm{argmax}\; n_l\f$ (index of the block with maximal \f$n_l\f$) */
  const Structure *getMaxBlock() const { 
    return getBlock(myMaxNlInd); 
  }
  /**@}*/
  
  bool isSameGamma() const { 
    return myIsSameGamma; 
  }
};


