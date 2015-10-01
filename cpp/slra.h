/*********************************
 * Header file fo SLRA package
 *********************************/

/* slra.h: SLRA header file */
#ifndef _SLRA_H_
#define _SLRA_H_

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_uint.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit_nlin.h> /* Levenberge-Marquardt */
#include <gsl/gsl_multimin.h>      /* BFGS Newton-type     */
}

#include "Exception.h"
#include "Cholesky.h"
#include "DGamma.h"
#include "Structure.h"
#include "MuDependentStructure.h"
#include "StationaryStructure.h"
#include "StripedStructure.h"
#include "StripedCholesky.h"
#include "StripedDGamma.h"
#include "HLayeredBlWStructure.h"
#include "HLayeredElWStructure.h"
#include "MuDependentCholesky.h"
#include "StationaryCholesky.h"
#include "StationaryCholeskySlicot.h"
#include "MuDependentDGamma.h"
#include "StationaryDGamma.h"
#include "PhiStructure.h"

#include "VarproFunction.h"
#include "NLSFunction.h"
#include "NLSVarpro.h"
#include "NLSVarproPsiXI.h"
#include "NLSVarproVecR.h"
#include "NLSVarproPsiVecR.h"

#include "OptimizationOptions.h"

#include "slra_common.h"
#include "Log.h"
#include "slralapack.h"

#include "MyIterationLogger.h"
#include "SLRAObject.h"


/** @defgroup MainFunctions
 * Global functions. */

/** \mainpage SLRA package C++ core library (efficient implementation of \cite slra-software)
 *
 * \section intro_sec Introduction
 *
 * This library provides efficient C++ implementation of the variable projection methods 
 * for structured low-rank approximation, described in \cite slra-efficient.
 *
 * The library is written in C++ in the object-oriented paradigm. The main mathematical objects
 * in \cite slra-efficient are represented as C++ objects. All the object methods and functions
 * are documented with references to correspinding formulas in \cite slra-efficient. The library
 * uses GSL, BLAS and LAPACK, and some other libraries.
 *
 * The library is a part of the SLRA package described in \cite slra-software. Thus the library
 *  provides interfaces to MATLAB (via mex-object) and R (via R package) languages.
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
 * are performed by VarproFunction object, which is constructed 
 * given Structure, \f$\Phi\f$ matrix, and rank reduction.
 *
 * VarproFunction object implements high-level algorithms for computing
 * cost function, gradient, pseudo-Jacobian and Jacobian given a general affine
 * structure.
 *
 * VarproFunction object also contains functions for computing initial approximation.
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



