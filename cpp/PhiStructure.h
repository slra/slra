/** Implementation Structure class for the structures of the form.
 * \f$\Phi \mathscr{S}'(p)\f$, there \f$\Phi\f$ is a full row rank matrix,
 * and \f$\mathscr{S}'(p)\f$ is a MuDependentStructure or StationaryStructure.
 */
class PhiStructure : public Structure {
private:
  Structure *myPStruct;
  gsl_matrix *myPhiT;
  gsl_matrix *myTempStMat;

  gsl_matrix *createPhiTRt( const gsl_matrix *Rt ) const;
public:
  /** Clolesky class for PhiStructure */
  class PhiCholesky : public Cholesky {
    friend class PhiStructure;
    const PhiStructure *myStruct;
    Cholesky *myParent;
    
  protected:
    /** Constructs a Cholesky object for PhiStructure.
     * @param[in] s    Pointer to the corresponding PhiStructure.
     * @param[in] d    number of rows \f$d\f$ of  the matrix \f$R\f$ */
    PhiCholesky( const PhiStructure *s, size_t d );
  public:
    virtual ~PhiCholesky();
    virtual void calcGammaCholesky( const gsl_matrix *Rt, double reg );
    virtual void multInvCholeskyVector( gsl_vector * yr, long trans );
    virtual void multInvGammaVector( gsl_vector * yr );
  };
  
  class PhiDGamma : virtual public DGamma {
    friend class PhiStructure;
    const PhiStructure *myStruct;
    DGamma *myParent;
    gsl_matrix *myPhi;
  protected:
    /** Constructs a Cholesky object for PhiStructure.
     * @param[in] s    Pointer to the corresponding PhiStructure.
     * @param[in] d    number of rows \f$d\f$ of  the matrix \f$R\f$ */
    PhiDGamma( const PhiStructure *s, size_t d );
  public:
    virtual ~PhiDGamma();
    virtual void calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt,
                               const gsl_matrix *Yt );
    virtual void calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt,
                                size_t j_1, size_t i_1, const gsl_vector *y,
                                const gsl_matrix *Phi = NULL );
  };


  /** Constructs PhiStructure from a pointer to a structure object.
   * @param PhiT the matrix \f$\Phi^{\top}\f$,
   *             where \f$\Phi \in \mathbb{R}^{m\times m'}\f$ is full row rank
   *             and \f$m'\f$ is the number of rows of the underlying structure.
   * @param S    a pointer to a  Structure object
   */
  PhiStructure( const gsl_matrix *PhiT, Structure *S );
  virtual ~PhiStructure();
  
  /** @name Implementing Structure interface */
  /**@{*/
  virtual size_t getN() const { return myPStruct->getN(); }
  virtual size_t getNp() const { return myPStruct->getNp(); }
  virtual size_t getM() const { return myPhiT->size2; }
  virtual void multByWInv( gsl_vector* p, long deg = 2 ) const {
    myPStruct->multByWInv(p, deg);
  }

  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ) ;
 
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt,
                                  const gsl_vector *y,
                                  double alpha = -1, double beta = 1,
                                  bool skipFixedBlocks = true );
  
  virtual Cholesky *createCholesky( size_t d ) const {
    return new PhiCholesky(this, d);
  }
  virtual DGamma *createDGamma( size_t d ) const {
    return new PhiDGamma(this, d);
  }
  /**@}*/
  
  /** @name StripedStructure-specific methods */
  /**@{*/
  /** Returns \f$\Phi^{\top}\f$ */
  const gsl_matrix *getPhiT() const { return myPhiT; }
  /**@}*/
};
