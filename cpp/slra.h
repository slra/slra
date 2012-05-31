/*********************************
 * Header file fo SLRA package
 *********************************/

/* slra.h: SLRA header file */
#ifndef _SLRA_H_
#define _SLRA_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h> /* Levenberge-Marquardt */
#include <gsl/gsl_multimin.h>      /* BFGS Newton-type     */
#include <gsl/gsl_blas.h>

#include "slra_optimization.h"
#include "slra_basic.h"
#include "slra_striped.h"

#include "slra_layered_hankel.h"
#include "slra_layered_hankel_weighted.h"
#include "slra_cholesky_btbanded.h"
#include "slra_dgamma_btbanded.h"

#include "slra_computation.h"

#include "slra_common.h"
#include "slralapack.h"


/** Function that creates appropriate Mosaic structure from weights and checks its
 *    consistence with \f$n_p\f$
 * @ingroup MainFunctions 
 * @param [in]     ml      Vector of \f$\bf m\f$ \sa MosaicHStructure
 * @param [in]     nk      Vector of \f$\bf n\f$ \sa MosaicHStructure  
 * @param [in]     wk      Vector of weights. If NULL, or of sizes \f$q\f$, \f$qN\f$ - 
 *                         MosaicHStructure is constructed, otherwise - WMosaicHStructure
 * @param [in]     d       Rank reduction  (\f$m-r\f$) 
 * @param [in]     np_comp Value of \f$n_p\f$
 */                
Structure *createMosaicStructure( gsl_vector * ml,  gsl_vector *nk, 
               gsl_vector * wk, int np_comp );



/** Main function that runs SLRA optimization
 * @ingroup MainFunctions 
 * @param [in]     p_in  Input parameter vector \f$p\f$  
 * @param [in]     s     Structure specification  
 * @param [in]     d     Rank reduction  
 * @param [in,out] opt   OptimizationOptions object
 * @param [in]     Rini  Matrix for initial approximation  
 * @param [in]     Phi   \f$\Phi\f$ matrix
 *                    (identity if <tt>Phi == NULL</tt> )   
 * @param [in]     Psi   \f$\Psi\f$ matrix
 *                       (identity if <tt>Psi == NULL</tt> )   
 * @param [out]    p_out Approximation \f$\widehat{p}\f$  
 *                       (not computed if <tt>p_out == NULL</tt> )
 * @param [out]    R_out Output parameter vector 
 *                       (not computed if <tt>R_out == NULL</tt> )
 * @param [out]    v_out Covariance matrix for X
 */
void slra( const gsl_vector *p_in, Structure* s, int d, 
          OptimizationOptions* opt, gsl_matrix *Rini, gsl_matrix *Phi, 
          gsl_matrix *Psi, gsl_vector *p_out, gsl_matrix *r_out, 
          gsl_matrix *v_out );

/** Main function that runs GSL optimization
 * @ingroup MainFunctions 
 * @param [in]     F     CostFunction object
 * @param [in,out] opt   Optimization options
 * @param [in,out] x_vec Vector containing initial approximation and returning
 *                       the minimum point 
 * @param [out]    v     Covariance matrix for x
 */
int gsl_optimize( CostFunction *F, OptimizationOptions *opt, 
                       gsl_vector* x_vec, gsl_matrix *v );

/** @defgroup MainFunctions
 * Global functions. */

/** \mainpage Using the SLRA C++ library
 *
 * \section intro_sec Introduction
 *
 * This library provides efficient implementation of the \c slra function
 * described in
 *
 * \verbatim
TECHREPORT{MarkovskyUsevich12-Software,
  author = {I. Markovsky and K. Usevich},
  title = {Software for weighted structured low-rank approximation},
  institution = {ECS, Univ. of Southampton},
  year = {2012},
  address = {\url{http://eprints.soton.ac.uk/}},
  texfile = {slra-ext.nw}
}
\endverbatim
 *
 * The library uses object-oriented paradigm. Most part 
 * of mathematical objects described in the manual for SLRA package
 * have the C++ counterpart objects. 
 *
 * The library makes extensive use GSL, BLAS and LAPACK libraries, so these
 * libraries should be linked when using the \c slra library.
 *
 * \section install_sec Entry point for the library
 * The main function for the library is slra(). The parameters for the
 * function are decribed below.
 * \copydetails slra()
 *
 * \subsection notes_sec Notes
 * \li All vectors and matrices are represented by (\c gsl_vector \c *) and 
 *     (\c gsl_matrix \c *) types from the GSL library. 
 * \li Important!  Note that transposed versions of matrices are used.
 *     This is caused by row-major storage of GSL library (compared to 
 *     standard column-major order in MATLAB and R).
 *
 * \subsection struct_sec  Structure specification 
 * Weights and structure are specified by creating a specific
 *    Structure object, which can be of following types:
 * \li LayeredHStructure - class for layered Hankel structure with blockwise weights
 * \li WLayeredHStructure - class for layered Hankel structure with elementwise weights
 * \li MosaicHStructure and WMosaicHStructure --- derivative classes for Mosaic Hankel
 *     structure, that are based on layered Hankel structure. 
 *
 * Function createMosaicStructure() helps to construct appropriate Mosaic Hankel structure 
 *     based on vectors  \f$\bf m\f$, \f$\bf n\f$ and a vector of weights.
 *
 *  Structure class represents and abstract class for general affine structure,
 *   and has the following methods to implement: determine sizes, fill in matrix, 
 *   compute correction for given \f$R\f$.
 *
 * \subsection cost_sec Cost function evaluation
 * Structure class can instantiate Cholesky and DGamma objects, which
 * represent respectively operations with Cholesky decomposition and derivatives
 * of the \f$\Gamma\f$ matrix for a given instance of structure. 
 * The Cholesky and DGamma instances are optimized for a specific structures.
 *
 * Computation of cost functions and derivatives
 * are performed by CostFunction object, which is constructed 
 * given Structure, \f$\Phi\f$ matrix, and rank reduction.
 *
 * CostFunction object implements high-level algorithms for computing
 * cost function, gradient, pseudo-Jacobian and Jacobian given a general affine
 * structure.
 *
 * CostFunction object also contains functions for computing initial approximation.
 *
 * \subsection opt_sec Optimization options
 * Optimization options are implemented in OptimizationOptions class,
 * and reflect all the options described in the manual. OptimizationOptions
 * documentation contains description of default values and helper function
 * for parsing \c opt.disp and \c opt.method option from a string.
 *
 * Optimization inside slra() is performed by gsl_optimize() function.
 */

#endif /* _SLRA_H_ */



