#ifndef _slra_mex_utils_h
#define _slra_mex_utils_h

#include <limits>
using namespace std;
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include "slra.h"

#ifndef BUILD_MEX_OCTAVE
#include "matrix.h"
#endif
#include "mex.h"

gsl_matrix M2trmat(  const mxArray * mat );

gsl_vector M2vec( const mxArray * mat );

char *M2Str( mxArray *myMat, char *str, size_t max_len );

#define MATStoreOption(MAT, opt, name, lvalue, uvalue)  \
  do {                                                   \
    mxArray *fld = mxGetField(MAT, 0, #name);               \
    if (fld != NULL && mxGetN(fld) != 0 && mxGetM(fld) != 0) { \
      double val = mxGetScalar(fld);	\
      if (val < lvalue) {					\
        mexWarnMsgTxt("Ignoring optimization option '"#name"' "	\
                      "because '"#name"' < "#lvalue".");	\
      } else if (val > uvalue) {			        \
        mexWarnMsgTxt("Ignoring optimization option '"#name"' "	\
                      "because '"#name"' > "#uvalue".");	\
      } else {                                                  \
        opt.name = val;                                         \
      }                                                         \
    } 								\
  } while (0)

#define STR_MAX_LEN 200

void mexFillOpt( const mxArray *Mopt, OptimizationOptions &opt, 
                 gsl_matrix & Rini, gsl_matrix &Psi, size_t m, size_t r );

#ifdef __cplusplus
extern "C" {
#endif
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
#ifdef __cplusplus
}
#endif

#endif /* _slra_mex_utils_h */