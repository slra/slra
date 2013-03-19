/** Prototype class for striped structure.
 * Abstract class that facilitates creation of structures of form
 * \f$\mathcal{S}(p) = 
 * \begin{bmatrix} 
 *   \mathcal{S}_1 (p^{(1)})\\ \vdots \\\mathcal{S}_N (p^{(N)})
 * \end{bmatrix}\f$
 * All structures must have the same \f$m\f$.
 */
class StripedStructure : public Structure {
  size_t myBlocksN;
  Structure **myStripe;

  /* Helper variables */
  size_t myN;
  size_t myNp;
  size_t myMaxNkInd;
  bool myIsSameGamma;
public:
  /** Constructs striped structure from array.
   * @param blocksN \f$N\f$ 
   * @param stripe array of Structure objects
   */
  StripedStructure( size_t blocksN, Structure *stripe[], 
                    bool isSameGamma = false );
  virtual ~StripedStructure();

  /** @name Implementing Structure interface */
  /**@{*/
  virtual size_t getN() const { return myN; }
  virtual size_t getNp() const { return myNp; }
  virtual size_t getM() const { return myStripe[0]->getM(); }
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p,
                                bool premultInvW = false ) ;
  virtual void correctP( gsl_vector* p, const gsl_matrix *R, 
                         const gsl_vector *yr, long wdeg = 2 );
  virtual Cholesky *createCholesky( size_t D ) const;
  virtual DGamma *createDGamma( size_t D ) const;
  /**@}*/
  
  /** @name Structure-specific methods */
  /**@{*/
  size_t getBlocksN() const { return myBlocksN; } /**< Returns \f$N\f$ */
  /** Returns \f$\mathcal{S}_{k}\f$ */
  const Structure *getBlock( size_t k ) const { 
    return myStripe[k]; 
  }
  /** Returns \f$\mathcal{S}_{k_{max}}\f$, where 
   * \f$k_{max} = argmax \{n_k\}_{k=1}^{N}\f$  */
  const Structure *getMaxBlock() const { 
    return getBlock(myMaxNkInd); 
  }
  /**@}*/
  
  bool isSameGamma() const { 
    return myIsSameGamma; 
  }
};


