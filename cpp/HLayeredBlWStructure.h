/** Layered Hankel structure with blockwise weights.
 * The layered Hankel structure is a structure of form 
 * \f$ \mathcal{H}_{{\bf m}, n} := 
 * \begin{bmatrix} \mathcal{H}_{n,m_1} & \cdots & \mathcal{H}_{n,m_q} 
 * \end{bmatrix} \f$  with \f$q\f$ blocks.
 */
class HLayeredBlWStructure : public StationaryStructure {
  typedef struct {
    size_t blocks_in_row;       /* Number of blocks in a row of Ci */
    double inv_w;            /* Square root of inverse of the weight */
  } Layer;

  size_t myQ;	                /* number of layers */
  size_t myN;
  size_t myM;
  size_t myMaxLag;
  gsl_matrix **myA;
  Layer *mySA;	/* q-element array describing C1,...,Cq; */  

  void computeStats();
  void computeWkParams(); 
protected:
  size_t nvGetNp() const { return (myN - 1) * myQ + myM; }  
public:
  /** Constructs layered Hankel structure.
   * @param m_l \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}\f$
   * @param q   \f$q\f$
   * @param w_l vector of weights 
   * \f${\bf w} =\begin{bmatrix} w_1 & \cdots & w_q \end{bmatrix}\f$.
   * If w_l == NULL then \f${\bf w}\f$ is set to be
   * \f$\begin{bmatrix}1&\cdots&1\end{bmatrix}\f$. 
   */
  HLayeredBlWStructure( const double *m_l, size_t q, size_t n, 
                     const double *w_l = NULL );
  virtual ~HLayeredBlWStructure();
  /** @name Implementing Structure interface */
  /**@{*/
  virtual size_t getM() const { return myM; }
  virtual size_t getN() const { return myN; }
  virtual size_t getNp() const { return nvGetNp(); }
  virtual Cholesky *createCholesky( size_t D ) const;
  virtual DGamma *createDGamma( size_t D ) const;
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ); 
  virtual void multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt, 
                                   const gsl_vector *y, 
                                   double alpha = -1, double beta = 1,
                                   bool skipFixedBlocks = true ); 
  virtual void multByWInv( gsl_vector* p, long deg = 2 );
  /**@}*/
 
  /** @name Implementing StationaryStructure interface */
  /**@{*/
  virtual size_t getS() const { return myMaxLag; }
  virtual void WkB( gsl_matrix *res, long k, const gsl_matrix *B ) const;
  virtual void AtWkB( gsl_matrix *res, long k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWkB, double beta = 0 ) const;
  virtual void AtWkV( gsl_vector *res, long k,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWkV, double beta = 0 ) const;
  /**@}*/
 
  /** @name Layered-Hankel-specific methods */
  /**@{*/
  /** Returns a pointer to Wk matrix.  @todo remove.*/
  const gsl_matrix *getWk( size_t l ) const { return myA[l]; } 
  /** Returns \f$q\f$. */
  size_t getQ() const { return myQ; }
  /** Returns \f$\max \{m_l\}_{l=1}^{q}\f$ */
  size_t getMaxLag() const { return myMaxLag; }
  /** Returns \f$m_l\f$ */
  size_t getLayerLag( size_t l ) const { return mySA[l].blocks_in_row; }
  /** Checks whether \f$w_l=\infty\f$ (block is fixed) */
  bool isLayerExact( size_t l ) const { return (mySA[l].inv_w == 0.0); }
  /** Returns \f$w_l^{-1}\f$ */
  double getLayerInvWeight( size_t l ) const { return mySA[l].inv_w; }
  /** Returns number of parameters in each block: \f$m_l +n- 1\f$ */
  size_t getLayerNp( size_t l ) const { return getLayerLag(l) + getN() - 1; }
  /**@}*/
};


