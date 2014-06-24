/** Implementation of Cholesky class for the MuDependentStructure.
 */
class MuDependentCholesky : public Cholesky {
protected:
  const MuDependentStructure *myStruct; /// Pointer to the Structure object
  size_t myD;                           /// \f$d\f$  
    
  size_t myMu_1;                   /// \f$\mu-1\f$
  size_t myDMu;                    /// \f$d\mu\f$  
  size_t myDN;                     /// \f$dn\f$  
  size_t myDMu_1;                  /// \f$d(\mu-1)\f$    
  gsl_matrix *myTempVijtRt;        /// Temporary storage for \f$\mathrm{V}_{\#ij} R^{\top}\f$
  gsl_matrix *myTempGammaij;       /// Temporary storage for \f$\Gamma_{\#ij}\f$ 

  /** The packed representation for $\f$\mathrm{L}_{\Gamma}\f$ 
   * and upper block-triangular part of \f$\Gamma(R)\f$.
   * Stored in column-major order. */
  double *myPackedCholesky;
protected:  

  /** The function computes the upper triangular part of \f$\Gamma(R)\f$
   * and puts it in MuDependentCholesky::myPackedCholesky.
   * The method has a faster implementation in StationaryCholesky.  
   *
   * @todo Consider rewriting this function in style of
   * StationaryCholesky::computeGammaUpperPart. Now this function is now too complicated 
   * due to mixing banded matrix representation in dpbtrf and row-major order in GSL. 
   * Although the function was extensively tested, it may be difficult to modify it.
   * 
   * @param[in] Rt        the matrix \f$R^{\top} \in \mathbb{R}^{m \times d}\f$.
   * @param[in] reg_gamma a regularization parameter \f$\gamma\f$. 
   */
  virtual void computeGammaUpperTrg( const gsl_matrix *Rt, double reg = 0 );
public:
  /** Constructs the MuDependentCholesky object.
   * @param[in] s    Pointer to the corresponding MuDependentStructure.
   * @param[in] d     number of rows \f$d\f$ of  the matrix \f$R\f$ */
  MuDependentCholesky( const MuDependentStructure *s, size_t d );
  virtual ~MuDependentCholesky();

  /** @name Implementing Cholesky interface */
  /**@{*/
  virtual void calcGammaCholesky( const gsl_matrix *Rt, double reg = 0 );
  virtual void multInvCholeskyVector( gsl_vector * y_r, long trans );
  virtual void multInvGammaVector( gsl_vector * y_r );
  /**@}*/

  /** @name Wrappers for MuDependentStructure methods */
  /**@{*/
  /** @copybrief Structure::getN() 
   * See Structure::getN().  */   
  size_t getN() const { return myStruct->getN(); }
  /** @copybrief Structure::getM() 
   * See Structure::getM().*/   
  size_t getM() const { return myStruct->getM(); }
  /** @copybrief MuDependentStructure::getMu()
   * See MuDependentStructure::getM(). */   
  size_t getMu() const { return myStruct->getMu(); }
  /**@}*/

  /** @name MuDependentCholesky-specific methods */
  /**@{*/
  /** Returns \f$d\f$ (the number of rows of \f$R\f$). */
  virtual size_t getD() const { return myD; }
  /**@}*/
};








