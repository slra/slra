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


void myMexErrorH( const char *reason, const char *F, int ln, int gsl_err );

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


class SLRAObject {
  Structure *myS;
  VarproFunction *myF; 
  gsl_error_handler_t *old_gsl_err_h;
  static size_t myObjCnt;
public:
  SLRAObject( gsl_vector p_in, gsl_vector ml, gsl_vector nk,
              gsl_matrix perm, gsl_vector wk, gsl_vector rvec, bool isgcd = false );
  virtual ~SLRAObject();
  
  Structure *getS() { return myS; }
  VarproFunction *getF() { return myF; }
};

#ifdef __cplusplus
extern "C" {
#endif
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] );
#ifdef __cplusplus
}
#endif