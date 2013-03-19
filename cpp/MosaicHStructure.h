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
};
