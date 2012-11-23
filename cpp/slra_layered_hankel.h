/** Layered Hankel structure with blockwise weights.
 * The layered Hankel structure is a structure of form 
 * \f$ \mathcal{H}_{{\bf m}, n} := 
 * \begin{bmatrix} \mathcal{H}_{n,m_1} & \cdots & \mathcal{H}_{n,m_q} 
 * \end{bmatrix} \f$  with \f$q\f$ blocks.
 */
class LayeredHStructure : public StationaryStructure {
  typedef struct {
    size_t blocks_in_row;       /* Number of blocks in a row of Ci */
    double inv_w;            /* Square root of inverse of the weight */
  } Layer;

  int myQ;	                /* number of layers */
  size_t myN;
  size_t myM;
  size_t myMaxLag;
  gsl_matrix **myA;
  Layer *mySA;	/* q-element array describing C1,...,Cq; */  

  void computeStats();
  void computeWkParams(); 
protected:
  int nvGetNp() const { return (myN - 1) * myQ + myM; }  
public:
  /** Constructs layered Hankel structure.
   * @param m_l \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}\f$
   * @param q   \f$q\f$
   * @param w_l vector of weights 
   * \f${\bf w} =\begin{bmatrix} w_1 & \cdots & w_q \end{bmatrix}\f$.
   * If w_l == NULL then \f${\bf w}\f$ is set to be
   * \f$\begin{bmatrix}1&\cdots&1\end{bmatrix}\f$. 
   */
  LayeredHStructure( const double *m_l, size_t q, int n, 
                     const double *w_l = NULL );
  virtual ~LayeredHStructure();
  /** @name Implementing Structure interface */
  /**@{*/
  virtual int getM() const { return myM; }
  virtual int getN() const { return myN; }
  virtual int getNp() const { return nvGetNp(); }
  virtual Cholesky *createCholesky( int D ) const;
  virtual DGamma *createDGamma( int D ) const;
  virtual void fillMatrixFromP( gsl_matrix* c, const gsl_vector* p ); 
  virtual void correctP( gsl_vector* p, gsl_matrix *R, gsl_vector *yr,
                         int wdeg = 2 );
  /**@}*/
 
  /** @name Implementing StationaryStructure interface */
  /**@{*/
  virtual int getS() const { return myMaxLag; }
  virtual void WkB( gsl_matrix *res, int k, const gsl_matrix *B ) const;
  virtual void AtWkB( gsl_matrix *res, int k, 
                      const gsl_matrix *A, const gsl_matrix *B, 
                      gsl_matrix *tmpWkB, double beta = 0 ) const;
  virtual void AtWkV( gsl_vector *res, int k,
                      const gsl_matrix *A, const gsl_vector *V, 
                      gsl_vector *tmpWkV, double beta = 0 ) const;
  /**@}*/
 
  /** @name Layered-Hankel-specific methods */
  /**@{*/
  /** Returns a pointer to Wk matrix.  @todo remove.*/
  const gsl_matrix *getWk( int l ) const { return myA[l]; } 
  /** Returns \f$q\f$. */
  int getQ() const { return myQ; }
  /** Returns \f$\max \{m_l\}_{l=1}^{q}\f$ */
  int getMaxLag() const { return myMaxLag; }
  /** Returns \f$m_l\f$ */
  int getLayerLag( int l ) const { return mySA[l].blocks_in_row; }
  /** Checks whether \f$w_l=\infty\f$ (block is fixed) */
  bool isLayerExact( int l ) const { return (mySA[l].inv_w == 0.0); }
  /** Returns \f$w_l^{-1}\f$ */
  double getLayerInvWeight( int l ) const { return mySA[l].inv_w; }
  /** Returns number of parameters in each block: \f$m_l +n- 1\f$ */
  int getLayerNp( int l ) const { return getLayerLag(l) + getN() - 1; }
  /**@}*/
};


/** Mosaic Hankel structure with blockwise weights.
 * The mosaic Hankel structure is a structure of form 
 * \f$ \mathcal{H}_{{\bf m}, {\bf n}} := 
 * \begin{bmatrix} 
 * \mathcal{H}_{n_1,m_1} & \cdots & \mathcal{H}_{n_1,m_q} \\
 * \vdots                &        & \vdots              \\
 * \mathcal{H}_{n_N,m_1} & \cdots & \mathcal{H}_{n_N,m_q} \\
 * \end{bmatrix} \f$  with \f$N q\f$ blocks.
 */
class MosaicHStructure : public StripedStructure {
  bool myWkIsCol;
protected:
  static Structure **allocStripe( gsl_vector *m_l, gsl_vector *n_k,  
                                  gsl_vector *w );
public:
  /** Constructs MosaicHStructure object 
   * by combining StripedStructure and LayeredHStructure.
   * @param m_l vector 
   * \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}\f$
   * @param n_k vector 
   * \f${\bf n} = \begin{bmatrix}n_1 & \cdots & n_N\end{bmatrix}\f$
   * @param w vector of weights \f${\bf w}\f$, defining
   * the matrix of blockwise weights 
   * \f$W = \begin{bmatrix} 
   * w_{1,1} & \cdots & w_{1,q} \\
   * \vdots  &        & \vdots  \\
   * w_{N,1} & \cdots & w_{N,q} \\
   * \end{bmatrix}\f$ in a following way
   * \f[{\rm vec}\ W^{\rm T} = 
   * \begin{cases}
   *   \begin{bmatrix}1&\cdots&1\end{bmatrix}, &  {\bf w} = NULL, \\
   *   \begin{bmatrix}1&\cdots&1\end{bmatrix} \otimes {\bf w}, &
   *   {\rm length}\ ({\bf w}) = q, \\
   *   {\bf w}, & \mbox{otherwise}. \\
   * \end{cases}\f]
   */
  MosaicHStructure( gsl_vector *m_l, gsl_vector *n_k, gsl_vector *w );
  virtual ~MosaicHStructure() {}
  /** @name Implementing Structure interface */
  /**@{*/
  virtual Cholesky *createCholesky( int r ) const;
  /**@}*/
};


