/** Abstract class for computations with derivarives of \f$\Gamma(R)\f$. */
class DGamma {
public:  
  virtual ~DGamma() {}
  /** Calculate the nonlinear part of the matrix gradient.
   * Calculates the matrix
   * \f$A = 2 \sum\limits_{i,j}^{n,n} Y_{:,j} Y_{:,i}^\top R V_{\#ij}\f$.
   *
   * This exactly the matrix \f$A\in\mathbb{R}^{d\times m}\f$ that satisfies
   * \f$\mathrm{trace}(AH^{\top}) = y^{\top} d \Gamma(R, H) y \f$.
   * (See also eqn. \f$(\nabla_{d\times m})\f$ and \f$(df(R,H))\f$
   * in \cite slra-efficient.)
   *
   * @param[out] At  transposed matrix \f$A^{\top}\in\mathbb{R}^{m\times d}\f$
   * @param[in]  Rt  transposed matrix \f$R^{\top}\in\mathbb{R}^{m\times d}\f$
   * @param[in]  Yt  matrix \f$Y^{\top} \in \mathbb{R}^{n\times d}\f$,
   *                 see eqn. \f$(Y)\f$ in \cite slra-efficient. 
   * */
  virtual void calcYtDgammaY( gsl_matrix *At, const gsl_matrix *Rt, 
                   const gsl_matrix *Yt ) = 0;

  /** Calculates a part of Jacobian/pseudo Jacobian
   *
   * Calculates 
   * \f$z \leftarrow \left((z^{(1)}_j) \otimes e^{(i)} + z^{(2)}_{ij}\right)\f$,
   * where the vectors \f$z^{(\cdot)}_{\cdot}\f$ are defined as
   * \f[
   * z^{(1)}_{j} = \sum_{k=1}^{n} \mathop\mathrm{col} \left(
   * \Phi_{j,:} V_{\#1k} R^{\top}Y_{:,k}, \ldots,
   * \Phi_{j,:} V_{\#nk} R^{\top}Y_{:,k}
   * \right),
   * \f]
   * and
   * \f[
   * z^{(2)}_{ij} = \sum_{k=1}^{n} Y_{i,k} \mathop\mathrm{col} \left(
   *  R V_{\#1k} \Phi_{j,:}^{\top}, \ldots, R  V_{\#nk} \Phi_{j,:}^{\top}
   * \right),
   * \f]
   * where \f$\Phi \in \mathbb{R}^{m'' \times m}\f$ is a given matrix.
   *
   * Note: for \f$\Phi = I_m\f$ these vectors are exactly the ones defined in
   * \f$(z^{(1)}_j)\f$ and \f$(z^{(2)}_{ij})\f$ in \cite slra-efficient.
   *
   * @param[out] z    the result
   * @param[in]  j_1  \f$0\f$-based index \f$j_1\f$, such that 
   *                  \f$j=j_1+1\f$  and \f$0 \le j_1 <m''\f$.
   * @param[in]  i_1  \f$0\f$-based index \f$i_1\f$, such that 
   *                  \f$i=i_1+1\f$  and \f$0 \le i_1 <d \f$.   
   * @param[in]  Rt  transposed matrix \f$R^{\top}\in\mathbb{R}^{m\times d}\f$.
   * @param[in]  y   vector \f$y \in \mathbb{R}^{dn}\f$.
   * @param[in]  Phi matrix \f$\Phi \in \mathbb{R}^{m'' \times m}\f$.
   *             If `w_vec == NULL` then \f$\Phi = I_m\f$.
   */
  virtual void calcDijGammaYr( gsl_vector *z, const gsl_matrix *Rt, 
                   size_t j_1, size_t i_1, const gsl_vector *y,
                   const gsl_matrix *Phi = NULL ) = 0;
};






