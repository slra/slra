/* MEX function for calling stls.c

% Author: Ivan Markovsky, Last modified: November 2004.
%
% Reference: I. Markovsky and S. Van Huffel "High-performance numerical algorithms 
% and software for structured total least squares", Journal of Computational and 
% Applied Mathematics, 2005
*/

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

/* default constants for the exit condition */

/* field names for opt */
#define EPSABS_STR "epsabs"
#define EPSREL_STR "epsrel"
#define EPSGRAD_STR "epsgrad"
#define MAXITER_STR "maxiter"
#define REGGAMMA_STR "reggamma"
#define DISP_STR "disp"
#define STR_MAX_LEN 25

/* values for OPT.DISP_STR */
#define NOTIFY_STR "notify"
#define FINAL_STR "final"
#define ITER_STR "iter"
#define OFF_STR "off"

/* field names for s */
#define ARRAY_STR "a"
#define NUM_ROLES_STR "k"

/* field names for info */
#define FMIN_STR "fmin"
#define ITER_STR "iter"
#define TIME_STR "time"

void m_to_gsl_matrix( gsl_matrix*, double* );
void gsl_to_m_matrix( double*, gsl_matrix* ); 

void tolowerstr( char * str ) {
  char *c;
  for (c = str; *c != '\0'; c++) {
    *c = tolower(*c);
  }
}

#define IfCheckAndStoreFieldBoundL(name, lvalue)   \
          if (! strcmp(field_name, #name)) {  \
            opt.name = mxGetScalar(mxGetFieldByNumber(prhs[4], 0, l)); \
            if (opt.name < lvalue) { \
              opt.name = SLRA_DEF_##name; \
              mexWarnMsgTxt("Ignoring optimization option '"#name"' because '"#name"' < "#lvalue"."); \
            } \
          }
                
#define IfCheckAndStoreFieldBoundLU(name, lvalue, uvalue)   \
          if (! strcmp(field_name, #name)) {  \
            opt.name = mxGetScalar(mxGetFieldByNumber(prhs[4], 0, l)); \
            if (opt.name < lvalue || opt.name > uvalue) { \
              opt.name = SLRA_DEF_##name; \
              mexWarnMsgTxt("Ignoring optimization option '"#name"' because '"#name"' < "#lvalue" or '"#name"' > "#uvalue"."); \
            } \
          }

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{
  gsl_matrix *x = NULL, *v = NULL;
  gsl_vector *p = NULL;
  gsl_vector_view vec_p;
  data_struct s;
  opt_and_info opt;
  char str_buf[STR_MAX_LEN];
  char *str_disp[] = {"",  "notify", "final", "iter", "off" };
  
  char meth_codes[] = "lqn";
  char submeth_codes_lm[] = "ls";
  char submeth_codes_qn[] = "b2pf";
  char submeth_codes_nm[] = "b2pf";
  char *submeth_codes[] = { submeth_codes_lm, submeth_codes_qn, submeth_codes_nm };
  int submeth_codes_max[] = { sizeof(submeth_codes_lm) / sizeof(submeth_codes_lm[0]) -1, 
          sizeof(submeth_codes_qn) / sizeof(submeth_codes_qn[0]) -1, 
          sizeof(submeth_codes_nm) / sizeof(submeth_codes_nm[0]) -1 };
  
  
  
  double *s_matr;
  int s_matr_cols;
  int l, i; 
  int n_plus_d, np;

  size_t n,d;
  char err_msg[100];
  int has_x;
  int col_p_vector = 1;
  
  /* ---------- */
  /* Input data */
  /* ---------- */

  if (nrhs < 2) {
    mexErrMsgTxt("Error: at least two parameters (p,s) are needed.");
  }
  

  /* check p */
  if (mxGetN(prhs[0]) == 1) { 
    np = mxGetM(prhs[0]);
  } else {
    if (mxGetM(prhs[0]) != 1) {
      mexErrMsgTxt("Error: p is neither a column nor a row vector.");
    }
    col_p_vector = 0;
    np = mxGetN(prhs[0]);
  }

  /* structure description prhs[1] */
  const mxArray* strArray;
  if (mxIsStruct(prhs[1])) {
    mxArray* field;
    /* in this case prhs[1] should have fields NUM_ROLES_STR and ARRAY_STR */
    if ((field = mxGetField(prhs[1], 0, NUM_ROLES_STR)) == NULL) {
      mexErrMsgTxt("Error in the structure specification : field " NUM_ROLES_STR " undefined.");
    }
    s.k = (int) mxGetScalar(field);

    if ((strArray = mxGetField(prhs[1], 0, ARRAY_STR)) == NULL) {
      mexErrMsgTxt("Error in the structure specification : field " ARRAY_STR " undefined.");
    }
  } else {
    /* in this case k = 1, and prhs[1] is the array */
    s.k = 1;
    strArray = prhs[1];
  }
  if (mxGetN(strArray) < 1 || mxGetN(strArray) > 4) {
    mexErrMsgTxt("Error in the structure specification : size(s." ARRAY_STR ",2) < 1 or > 4.");
  }
  if (mxGetM(strArray) < 1 || mxGetM(strArray) > 10) {
    mexErrMsgTxt("Error in the structure specification : size(s." ARRAY_STR ",1) < 1 or > 10.");
  }
  s_matr = mxGetPr(strArray);
  s.q = mxGetM(strArray);
  s_matr_cols = mxGetN(strArray);

  /* Create s and check structure specification */
  n_plus_d = 0;
  for (l = 0; l < s.q; l++) {
    s.a[l].blocks_in_row = *(s_matr + l);
    s.a[l].nb = (s_matr_cols > 1) ? *(s_matr + s.q + l): 1;
    s.a[l].exact = (s_matr_cols > 2) ? *(s_matr + 2 * s.q + l): 0;
    s.a[l].toeplitz = (s_matr_cols > 3) ? *(s_matr + 3 * s.q + l): 0;
    n_plus_d += s.a[l].blocks_in_row * s.a[l].nb;
  }

  /* Get r (rank) */
  if (nrhs > 2) {
    n = (int)(mxGetScalar(prhs[2]));
    d = n_plus_d - n;
  } else {
    d = 1;
    n = n_plus_d - d;
  }
    
  /* Check initial approximation prhs[3] */
  has_x = (nrhs >= 4) && (mxGetM(prhs[3]) > 0);
  if (has_x) {
    /* check dimensions of prhs[3] */
    if ( n != mxGetM(prhs[3]) ) {/* check n */
      mexErrMsgTxt("Error: n ~= size(x,1).");
    }
    if ( d != mxGetN(prhs[3]) ) { /* check d */
      mexErrMsgTxt("Error: d ~= size(x,1).");
    }
  } 

  /* ---------------------- */
  /* Allocate and copy data */
  /* ---------------------- */

  /* convert x in GSL format and store it in xh */
  x = gsl_matrix_alloc(n, d);
  if (has_x) {
    m_to_gsl_matrix(x, mxGetPr( prhs[3] ) );
  }
  
  /* Allocate and copy parameters vector */
  p = gsl_vector_alloc(np);
  vec_p = gsl_vector_view_array(mxGetPr(prhs[0]), np);
  gsl_vector_memcpy(p, &vec_p.vector);

  v = gsl_matrix_alloc(n*d, n*d);
  
  /* optimization options prhs[4] */
  /* default options */
  slraAssignDefOptValues(opt);

  /* user supplied options */
  if (nrhs >= 5) {
    if (! mxIsStruct(prhs[4]))
      mexWarnMsgTxt("Ignoring 'opt'. The optimization options should be passed in a structure.");
    else {
      int nfields = mxGetNumberOfFields(prhs[4]);
      const char *field_name_ptr;
      char field_name[20], *c;
      for (l = 0; l < nfields; l++) {
	field_name_ptr = mxGetFieldNameByNumber(prhs[4], l);
	strcpy(field_name, field_name_ptr);
	tolowerstr(field_name);
	
	/* which option */
	if (! strcmp(field_name, DISP_STR)) {
 	  mxGetString(mxGetFieldByNumber(prhs[4], 0, l), str_buf, STR_MAX_LEN); 
	  tolowerstr(str_buf);
	    
	  for (i = sizeof(str_disp)/sizeof(str_disp[0]) - 1; i > 0; i--)  {
	    if (!strcmp(str_buf, str_disp[i])) {
	      opt.disp = i;
	      break;
	    }
	  }
	  if (i <= 0) {
	    mexWarnMsgTxt("Ignoring optimization option 'disp'. Unrecognized value.");
	  }
 	} else if (! strcmp(field_name, "method")) {
 	  mxGetString(mxGetFieldByNumber(prhs[4], 0, l), str_buf, STR_MAX_LEN); 
	  tolowerstr(str_buf);
	  
	  for (i = sizeof(meth_codes)/sizeof(meth_codes[0]) - 1; i >= 0; i--)  {
	    if (str_buf[0] == meth_codes[i]) {
	      opt.method = i;
	      break;
	    }
	  } 
	  
	  if (i < 0)  {
	    mexWarnMsgTxt("Ignoring optimization option 'method'. Unrecognized value.");
	    slraAssignDefOptValue(opt, method);
	    slraAssignDefOptValue(opt, submethod);
	  }


	  for (i = submeth_codes_max[opt.method] - 1; i >= 0; i--)  {
	    if (str_buf[1] == submeth_codes[opt.method][i]) {
	      opt.submethod = i;
	      break;
	    }
	  } 
	  if (i < 0)  {
	    mexWarnMsgTxt("Unrecognized or unspecified submethod - using default.");
	    slraAssignDefOptValue(opt, submethod);
	  }

 	} else IfCheckAndStoreFieldBoundL(maxiter, 0)  
        else IfCheckAndStoreFieldBoundLU(epsabs, 0, 1) 
        else IfCheckAndStoreFieldBoundLU(epsrel, 0, 1) 
        else IfCheckAndStoreFieldBoundLU(epsgrad, 0, 1) 
        else IfCheckAndStoreFieldBoundLU(epsx, 0, 1) 
        else IfCheckAndStoreFieldBoundLU(step, 0, 1) 
        else IfCheckAndStoreFieldBoundLU(tol, 0, 1) 
	else IfCheckAndStoreFieldBoundL(reggamma, 0) 
        else { 
 	  sprintf(err_msg, "Ignoring unrecognized optimization option '%s'.", field_name); 
	  mexWarnMsgTxt(err_msg); 
 	} 
      }
    }
  }
  
  /* --------------- */
  /* Call the solver */
  /* --------------- */

  slra(p, &s, x, v, &opt, has_x, (nlhs > 3));

  /* ------------------ */
  /* Assign the outputs */
  /* ------------------ */

  /* estimate plhs[0] */
  plhs[0] = mxCreateDoubleMatrix(n, d, mxREAL);
  gsl_to_m_matrix(mxGetPr( plhs[0] ), x);

  /* output info */
  if (nlhs > 1) {
    l = 1;
    const char *field_names[] = {FMIN_STR, ITER_STR, TIME_STR};
    mxArray *temp;
    plhs[1] = mxCreateStructArray(1, &l, 3, field_names);
    temp = mxCreateDoubleScalar( opt.fmin ); 
    mxSetField(plhs[1], 0, FMIN_STR, temp);
    temp = mxCreateDoubleScalar( opt.iter ); 
    mxSetField(plhs[1], 0, ITER_STR, temp);
    temp = mxCreateDoubleScalar( opt.time ); 
    mxSetField(plhs[1], 0, TIME_STR, temp);
  }

  /* covariance matrix */
  if (nlhs > 2) {
    plhs[2] = mxCreateDoubleMatrix(n*d, n*d, mxREAL);
    gsl_to_m_matrix(mxGetPr( plhs[2] ), v);
  }
  
  /* compute correction */
  if (nlhs > 3) {
    if (col_p_vector) {
      plhs[3] = mxCreateDoubleMatrix(np, 1, mxREAL);
    } else {
      plhs[3] = mxCreateDoubleMatrix(1, np, mxREAL);
    }
    vec_p = gsl_vector_view_array(mxGetPr(plhs[3]), np);
    gsl_vector_memcpy(&vec_p.vector, p);
  }


  /* --------------------- */
  /* Free allocated memory */
  /* --------------------- */

  gsl_vector_free(p);
  gsl_matrix_free(x);
  gsl_matrix_free(v);
}






