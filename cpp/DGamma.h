/** Abstract class for computations with derivarives of \f$\Gamma(R)\f$. */
class DGamma {
public:  
  virtual ~DGamma() {}
  /** Calculate the nonlinear part of the matrix gradient.
   * Computes the matrix \f$A\in\mathbb{R}^{d\times m}\f$ such that
   * \f$\mathrm{trace}(AH^{top}) = y^{\top} d \Gamma(R, H) y \f$.
   * See also eqn. \f$(df(R,H))\f$ in \cite slra-efficient.
   *
   * @param[out] Bt  transposed matrix \f$B^{\top}\in\mathbb{R}^{m\times d}\f$
   * @param[in]  Rt  transposed matrix \f$R^{\top}\in\mathbb{R}^{m\times d}\f$
   * @param[in]  Yt  matrix \f$Y^{\top} \in \mathbb{R}^{n\times d}\f$,
   *                 see eqn. \f$(Y)\f$ in \cite slra-efficient. 
   * */
  virtual void calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt, 
                   const gsl_matrix *Yt ) = 0;

  /** Calculates \f$z \leftarrow \frac{\partial\Gamma}{\partial R_{ij}} y\f$.
   * See also eqn. \f$(z_{ij})\f$ in \cite slra-efficient. 
   * @param[out] z    the result of the multiplcation
   * @param[in]  j_1  \f$0\f$-based index \f$j_1\f$, such that 
   *                  \f$j=j_1+1\f$  and \f$0 \le j_1 <m \f$.   
   * @param[in]  i_1  \f$0\f$-based index \f$i_1\f$, such that 
   *                  \f$i=i_1+1\f$  and \f$0 \le i_1 <d \f$.   
   * @param[in]  Rt  transposed matrix \f$R^{\top}\in\mathbb{R}^{m\times d}\f$
   * @param[in]  y   vector \f$y \in \mathbb{R}^{dn}\f$
   * */
  virtual void calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt, 
                   size_t j_1, size_t i_1, const gsl_vector *y ) = 0;
};






