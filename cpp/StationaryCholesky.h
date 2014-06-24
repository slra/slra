/** Implementation of Cholesky class for the StationaryStructure.
 * A descendant of MuDependentCholesky. Basically, only
 * MuDependentCholesky::computeGammaUpperTrg is reimplemented
 */
class StationaryCholesky : public MuDependentCholesky {
public:
  /** Constructs the StationaryCholesky object.
   * @param[in] s    Pointer to the corresponding StationaryStructure.
   * @param[in] d     number of rows \f$d\f$ of  the matrix \f$R\f$ */
  StationaryCholesky( const StationaryStructure *s, size_t d );
  virtual ~StationaryCholesky();

protected:
  const StationaryStructure *myStStruct;
  /** A temporary object for storing the first block row of \f$\Gamma(R)\f$, more precisely
   * \f[ \begin{bmatrix} \Gamma_0 & \Gamma_1 & \cdots & \Gamma_{\mu-1} & 0 \end{bmatrix},\f]
   * where \f$ \Gamma_k := R^{\top} \mathrm{V}_{k} R\f$.
   */
  gsl_matrix *myGammaK;   

  /** Reimplements  MuDependentCholesky::computeGammaUpperTrg() using the
   * fact that \f$\Gamma_{\#ij} = \Gamma_{j-i}\f$.
   *
   * @param[in] Rt        the matrix \f$R^{\top} \in \mathbb{R}^{m \times d}\f$.
   * @param[in] reg       a regularization parameter \f$\gamma\f$. 
   */
  virtual void computeGammaUpperTrg( const gsl_matrix *Rt, double reg = 0 );

  /** Computes all \f$\Gamma_k\f$ and puts them in StationaryCholesky::myGammaK.
   * @copydetails StationaryCholesky::computeGammaUpperTrg */
  virtual void computeGammak( const gsl_matrix *Rt, double reg = 0 );
};


