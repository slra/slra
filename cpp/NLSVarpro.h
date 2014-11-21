#ifndef _NLSVarpro_h
#define _NLSVarpro_h

/** A subclass of the NLSFunction class for VARPRO functions
 * Used for functions that come from application of VARPRO to SLRA.
 */
class NLSVarpro : public NLSFunction {
public:
    /** Returns \f$n\f$ --- the number of variables */
    virtual size_t getNvar() = 0;
    /** Returns \f$n_s\f$ --- the number of squares (dimension of \f$g\f$) */
    virtual size_t getNsq() = 0;
    
    /** Computes function \f$f\f$ and gradient */
    virtual void computeFuncAndGrad( const gsl_vector* x, double* f,
                                    gsl_vector *grad ) = 0;
    /** Computes the vector \f$g\f$ and the Jacobian (or pseudo-jacobian) */
    virtual void computeFuncAndJac( const gsl_vector* x, gsl_vector *res,
                                   gsl_matrix *jac ) = 0;
    
    virtual size_t getD() = 0;
    virtual size_t getM() = 0;
    
    virtual void RTheta2x( gsl_matrix *RTheta, gsl_vector *x ) = 0;
    virtual void x2RTheta( gsl_matrix *RTheta, const gsl_vector *x ) = 0;
    
    /* To remove: */
    virtual void computePhat( gsl_vector* p, const gsl_vector* x ) = 0;
    virtual void computeDefaultx( gsl_vector *x ) = 0;
    
    virtual gsl_matrix x2xmat( const gsl_vector *x ) = 0;
};

#endif
