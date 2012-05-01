/* MEX function for calling stls.c */
#include <limits>
using namespace std;
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>
#include "slra.h"

#ifndef BUILD_MEX_OCTAVE
#include "matrix.h"
#endif
#include "mex.h"

/* field names for opt */
#define DISP_STR "disp"
#define RANK_STR "r"
#define RINI_STR "Rini"
#define PERM_STR "phi"
#define WK_STR "w"

/* field names for s */
#define STR_ML "m"
#define STR_NK "n"

/* field names for info */
#define RH_STR "Rh"
#define VH_STR "Vh"
#define FMIN_STR "fmin"
#define ITER_STR "iter"
#define TIME_STR "time"
#define METHOD_STR "method"


void SLRA_mex_error_handler( const char * reason, const char * file, int line, 
         int gsl_errno) {
  char err_msg[250];
  
  throw new Exception("GSL error #%d at %s:%d:  %s", file, line, gsl_errno, 
                reason);
}

gsl_matrix M2trmat( mxArray * mat ) {
  gsl_matrix res = { 0, 0, 0, 0, 0, 0 };
  if (mat != NULL && mxGetN(mat) != 0 && mxGetM(mat) != 0) {
    res = gsl_matrix_view_array(mxGetPr(mat),mxGetN(mat),mxGetM(mat)).matrix;
  }
  return res;
}

gsl_vector M2vec( const mxArray * mat ) {
  gsl_vector res = { 0, 0, 0, 0, 0 };
  if (mat != NULL && mxGetN(mat) != 0 && mxGetM(mat) != 0) {
    res = gsl_vector_const_view_array(mxGetPr(mat), mxGetN(mat)*mxGetM(mat)).vector;
  }
  return res;
}

char *M2Str( mxArray *myMat, char *str, int max_len ) {
  if (myMat == NULL) {
    *str = 0;
  } else {   
    mxGetString(myMat, str, max_len); 
  }  
  return str;
}

#define MATStoreOption(MAT, opt, name, lvalue, uvalue)  \
  do {                                                   \
    mxArray *fld = mxGetField(MAT, 0, #name);               \
    if (fld != NULL) {			         	\
      opt.name = mxGetScalar(fld);	\
      if (opt.name < lvalue) {					\
        opt.name = SLRA_DEF_##name;				\
        mexWarnMsgTxt("Ignoring optimization option '"#name"' "	\
                      "because '"#name"' < "#lvalue".");	\
      } else if (opt.name > uvalue) {			        \
        opt.name = SLRA_DEF_##name;				\
        mexWarnMsgTxt("Ignoring optimization option '"#name"' "	\
                      "because '"#name"' > "#uvalue".");	\
      }                                                         \
    } else { \
      opt.method =  SLRA_DEF_##disp; \
    }								\
  } while (0)

#define STR_MAX_LEN 200

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  gsl_error_handler_t *new_gsl_err_h = SLRA_mex_error_handler;
  gsl_error_handler_t *old_gsl_err_h = gsl_set_error_handler(new_gsl_err_h);
  char str_buf[STR_MAX_LEN];
  gsl_matrix rini = { 0, 0, 0, 0, 0, 0 }, rh_view = { 0, 0, 0, 0, 0, 0 },
             vh_view = { 0, 0, 0, 0, 0, 0 };
  Structure *myStruct;
  int  m, r, was_error = 0;
  opt_and_info opt;
  
  /* Input data */
  if (nrhs < 3) {
    mexErrMsgTxt("Error: at least two parameters (p, s, r) are needed.");
  }

  try {
    const gsl_vector p_in = M2vec(prhs[0]);
    gsl_vector ml = M2vec(mxGetField(prhs[1], 0, STR_ML));
    gsl_vector nk = M2vec(mxGetField(prhs[1], 0, STR_NK));
    if (ml.size == 0 || nk.size == 0) {
      mexErrMsgTxt("s.m and s.n should be nonempty vectors");   
    }

    int np = compute_np(&ml, &nk);
    if (p_in.size < np) {
      mexErrMsgTxt("Size of vector p less that the structure requires");   
    } else if (p_in.size > np) {
      mexWarnMsgTxt("Size of vector p more that the structure requires");   
    } 
    gsl_matrix perm = M2trmat(mxGetField(prhs[1], 0, PERM_STR));
    gsl_vector wk = M2vec(mxGetField(prhs[1], 0, WK_STR));
    r = mxGetScalar(prhs[2]);

    myStruct = createMosaicStructure(&ml, &nk, vecCheckNULL(wk), p_in.size);
    m = perm.data == NULL ? myStruct->getNplusD() : perm.size2;
    if (r <= 0 || r >= m) {
      throw new Exception("Incorrect rank\n");   
    }

    /* user supplied options */
    AssignDefOptValues(opt);
    if (nrhs > 3) {
      if (! mxIsStruct(prhs[3])) {
        mexWarnMsgTxt("Ignoring 'opt'. The optimization options "
	              "should be passed in a structure.");
      } else {
        rini = M2trmat(mxGetField(prhs[3], 0, RINI_STR));
        if (rini.data != NULL && (rini.size2 != (m - r) || rini.size1 != m)) {
          throw new Exception("Incorrect Rini\n");   
        }
        opt.disp = String2Disp(M2Str(mxGetField(prhs[3], 0, DISP_STR), 
                                     str_buf, STR_MAX_LEN));
        String2Method(M2Str(mxGetField(prhs[3], 0, METHOD_STR), 
                            str_buf, STR_MAX_LEN), &opt);
        MATStoreOption(prhs[3], opt, maxiter, 0, numeric_limits<int>::max());
        MATStoreOption(prhs[3], opt, epsabs, 0, 1);
        MATStoreOption(prhs[3], opt, epsrel, 0, 1);
        MATStoreOption(prhs[3], opt, epsgrad, 0, 1);
        MATStoreOption(prhs[3], opt, epsx, 0, 1);
        MATStoreOption(prhs[3], opt, step, 0, 1);
        MATStoreOption(prhs[3], opt, tol, 0, 1);
        MATStoreOption(prhs[3], opt, reggamma, 0, 
                       numeric_limits<double>::max());
        MATStoreOption(prhs[3], opt, ls_correction, 0, 1);
        MATStoreOption(prhs[3], opt, gcd, 0, 1);
      }
    }  
  
    /* output info */
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    gsl_vector p_out = M2vec(plhs[0]);

    if (nlhs > 1) {
      int l = 1;
      const char *field_names[] = { RH_STR, VH_STR, FMIN_STR, ITER_STR, 
                                    TIME_STR };
      plhs[1] = mxCreateStructArray(1, &l, 5, field_names);
    
      mxArray *rh = mxCreateDoubleMatrix((m - r), m, mxREAL);
      mxArray *vh = mxCreateDoubleMatrix((m - r) * r, (m - r) * r, mxREAL);
      rh_view = M2trmat(rh);
      vh_view = M2trmat(vh);
      mxSetField(plhs[1], 0, RH_STR, rh);
      mxSetField(plhs[1], 0, VH_STR, vh);
    }

    slra(vecCheckNULL(p_in), myStruct, r, &opt, matCheckNULL(rini), matCheckNULL(perm), 
         vecCheckNULL(p_out), matCheckNULL(rh_view), matCheckNULL(vh_view));

    if (nlhs > 1) {
      mxSetField(plhs[1], 0, FMIN_STR, mxCreateDoubleScalar(opt.fmin));
      mxSetField(plhs[1], 0, ITER_STR, mxCreateDoubleScalar(opt.iter));
      mxSetField(plhs[1], 0, TIME_STR, mxCreateDoubleScalar(opt.time));
    }
  } catch (Exception *e) {
    strncpy(str_buf, e->getMessage(), STR_MAX_LEN - 1);
    str_buf[STR_MAX_LEN - 1] = 0;
    was_error = 1;
    delete e;
  } 

  gsl_set_error_handler(old_gsl_err_h);
  if (myStruct != NULL) {
    delete myStruct;
  }
  if (was_error) {
    mexErrMsgTxt(str_buf);
  }
}






