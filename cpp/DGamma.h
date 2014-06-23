/** Abstract class for differentiating \f$\Gamma(R)\f$.
 * Computations with \f$\Gamma(R)\f$ */
class DGamma {
public:  
  virtual ~DGamma() {}
  /** Calculate the nonlinear part of the gradient.
   * Updates the gradient \f$grad \leftarrow grad + A\f$, where
   * \f$A\f$ is defined by \f$trace(A,H) = y_r d \Gamma(R, H) y_r\f$
   * */
  virtual void calcYrtDgammaYr( gsl_matrix *grad, const gsl_matrix *R, 
                   const gsl_vector *yr ) = 0;
  /** Calculate the nonlinear part of pseudojacobian.
   * Calculates \f$\displaystyle res \leftarrow \frac{\partial}{\partial R_{ij}} 
    \Gamma \left(R \right)y_r\f$
   * */
  virtual void calcDijGammaYr( gsl_vector *res, const gsl_matrix *R, 
                   size_t i, size_t j, const gsl_vector *Yr ) = 0;
};






