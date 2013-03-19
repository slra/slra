/** Mosaic Hankel structure with elementwise weights.
 * @copydetails MosaicHStructure
 */
class WMosaicHStructure : public StripedStructure {
protected:
  static Structure **allocStripe( gsl_vector *m_l, gsl_vector *n_k,  
                gsl_vector *w );
public:
  /** Constructs WMosaicHStructure object.
   * @param m_l \f${\bf m} = \begin{bmatrix}m_1 & \cdots & m_q\end{bmatrix}\f$
   * @param n_k \f${\bf n} = \begin{bmatrix}n_1 & \cdots & n_N\end{bmatrix}\f$
   * @param w vector of weights 
   * \f${\bf w} =\begin{bmatrix} w_1 & \cdots & w_{n_p} \end{bmatrix}\f$.
   * If w_l == NULL then \f${\bf w}\f$ is set to be
   * \f$\begin{bmatrix}1&\cdots&1\end{bmatrix}\f$. 
   */
  WMosaicHStructure( gsl_vector *m_l, gsl_vector *n_k, gsl_vector *w );
  virtual ~WMosaicHStructure() {}
};


