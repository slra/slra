/* MEX function for calling stls.c

[ xh, info, vxh ] = stls( a, b, s, x0, opt )

Input arguments:

A, B - data.
S    - structure specification, a matrix with q rows and 3 columns:
  S(i,1) - structure type [ 1 | 2 | 3 | 4 ] of the i-th block,
    1 - Toeplitz, 2 - Hankel, 3 - Unstructured, 4 - Exact (noise free),
  S(i,2) - number of columns of the i-th block,
  S(i,3) - number of columns of the blocks in the i-th block.
For block-Toeplitz/Hankel strucutred problem with K x nu(i) size blocks, 
S is a structure with fiels K and A, where S.A is the q x 2 matrix described 
above and S.K = K.
X0   - initial approximation (default TLS).
OPT  - optimization parameters, structure with fields: 
  OPT.MAXITER - maximum number of iterations, 
  OPT.EPSREL, OPT.EPSABS, OPT.EPSGRAD - convergence tolerances, and 
  OPT.REGGAMMA - regularization parameter for gamma, absolute (to be changed)
  OPT.DISP - level of display ['notify','final','iter',off] (default 'notify').
Exit condition: # iterations >= MAXITER, |xi(k+1)-xi(k)| < EPSABS + EPSREL*|xi(k+1)|
for all i, where x(k) is the approximation on the i-th step, or ||f'|| < EPSGRAD, 
where f' is the gradient of the cost function.

Output arguments:

XH   - STLS estimate.
INFO - information on exit, structure with fields ITER, TIME, and FMIN:
  INFO.ITER - number of iterations for convergence,
  INFO.TIME - execution time,
  INFO.FMIN - value of the cost function.
VXH  - assymptotic covariance matrix of the estimate XH.

Note: The program can not treat the case length(P) < size(A,1) * size(B,2).

Reference: I. Markovsky and S. Van Huffel "High-performance numerical algorithms 
and software for structured total least squares", Journal of Computational and 
Applied Mathematics, 2005

Author: Ivan Markovsky, Last modified: November 2004.
*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>
#include "stls.h"

#ifndef MEX_OCTAVE
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
#define NUM_ROWS_STR "m"
#define RANK_REDUCTION_STR "d"

/* field names for info */
#define FMIN_STR "fmin"
#define ITER_STR "iter"
#define TIME_STR "time"

void m_to_gsl_matrix( gsl_matrix*, double* );
void gsl_to_m_matrix( double*, gsl_matrix* ); 



#define DEF_disp 3 
#define DEF_epsabs 0
#define DEF_epsrel 1e-6
#define DEF_epsgrad 1e-6
#define DEF_maxiter 100
#define DEF_reggamma 0.00001



#define IfCheckAndStoreFieldBoundL(name, lvalue)   \
          if (! strcmp(field_name, #name)) {  \
            opt.name = mxGetScalar(mxGetFieldByNumber(prhs[4], 0, l)); \
            if (opt.name < lvalue) { \
              opt.name = DEF_##name; \
              mexWarnMsgTxt("Ignoring optimization option '"#name"' because '"#name"' < "#lvalue"."); \
            } \
          }
                
#define IfCheckAndStoreFieldBoundLU(name, lvalue, uvalue)   \
          if (! strcmp(field_name, #name)) {  \
            opt.name = mxGetScalar(mxGetFieldByNumber(prhs[4], 0, l)); \
            if (opt.name < lvalue || opt.name > uvalue) { \
              opt.name = DEF_##name; \
              mexWarnMsgTxt("Ignoring optimization option '"#name"' because '"#name"' < "#lvalue" or '"#name"' > "#uvalue"."); \
            } \
          }




void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{
  gsl_matrix *a, *b, *x, *v;
  gsl_vector *p;
  gsl_vector_view vec_p;
  data_struct s;
  opt_and_info opt;
  char str_buf[STR_MAX_LEN];
  char str_codes[] = " THUE";
  char *str_disp[] = {"",  "notify", "final", "iter", "off" };
  

  double *s_m;
  int l, i; 
  int n_plus_d, np;

  size_t m, n, d;
  char err_msg[100];
  int has_ab, has_p, has_x;
  
  

  /* ---------- */
  /* Input data */
  /* ---------- */

  if (nrhs < 3) {
    mexErrMsgTxt("Error: at least three parameters (a,b,s) are needed.");
  }

  has_ab = (!mxIsEmpty(prhs[0])) &&  (!mxIsEmpty(prhs[1]));
  has_p = (nrhs >= 6) && (!mxIsEmpty(prhs[5]));
  has_x = nrhs >= 4 && !mxIsEmpty(prhs[3]);

  if (!(has_ab || has_p)) {
    mexErrMsgTxt("Either (a,b) or (p) should be given.");
  
  } 


  /* structure description prhs[2] */
  if (mxIsStruct(prhs[2])) {
    mxArray* field;
    /* in this case prhs[3] should have fields NUM_ROLES_STR and ARRAY_STR */
    if ((field = mxGetField(prhs[2], 0, NUM_ROLES_STR)) == NULL) {
      mexErrMsgTxt("Error in the structure specification : field " NUM_ROLES_STR " undefined.");
    }
    s.k = (int) mxGetScalar(field);

    if ((field = mxGetField(prhs[2], 0, NUM_ROWS_STR)) != NULL) {
      s.m = (int) mxGetScalar(field);
    } else {
      if (!has_ab) {
        mexErrMsgTxt("Error in the structure specification : number of rows (s.m) not specified.");
      }
    }

    if ((field = mxGetField(prhs[2], 0, RANK_REDUCTION_STR)) != NULL) {
      d = (int) mxGetScalar(field);
    } else {
      if (!( has_ab|| has_x)) {
        mexErrMsgTxt("Error in the structure specification : rank reduction (s.d) should specified if none of (a,b,x) is given.");
      }
    }

    
    if ((field = mxGetField(prhs[2], 0, ARRAY_STR)) == NULL) {
      mexErrMsgTxt("Error in the structure specification : field " ARRAY_STR " undefined.");
    }
    if (mxGetN(field) != 3) {
      mexErrMsgTxt("Error in the structure specification : size(s." ARRAY_STR ",2) ~= 3.");
    }
    s_m = mxGetPr( field );
    s.q = mxGetM( field );
  } else {
    /* in this case k = 1, and prhs[3] is the array */
    s.k = 1;
    if (mxGetN( prhs[2] ) != 3)
      mexErrMsgTxt("Error in the structure specification : size(s,2) ~= 3.");
    s_m = mxGetPr( prhs[2] );    
    s.q = mxGetM( prhs[2] );
  }


  n_plus_d = 0;
  /* creat s */
  for (l = 0; l < s.q; l++) {
    if (*(s_m+l) < 1 ||  *(s_m+l) > 4) {
      sprintf(err_msg, "Error: invalid structure specification '%d'.", (int) *(s_m+l));
      mexErrMsgTxt(err_msg);
    } else {
      s.a[l].type = str_codes[(int)(*(s_m+l))]; 
    }
    s.a[l].ncol = *(s_m + s.q + l);
    n_plus_d += s.a[l].ncol;
    s.a[l].nb   = *(s_m + 2*s.q + l);
  }


  if (!has_ab) { 
    m = s.m;
    if (has_x){
      n = mxGetM(prhs[3]);
      d = mxGetN(prhs[3]);
    } else {
      n = n_plus_d - d;
    }
    
   } else  {
    m = mxGetM( prhs[0] );
    s.m = m;
    n = mxGetN( prhs[0] );
    d = mxGetN( prhs[1] );

    /* check dimensions of a and b */
    if ( m != mxGetM(prhs[1]) )  { /* check m */
      mexErrMsgTxt("Error: size(a,1) ~= size(b,1).");
    }
  }

  if ( m < n ) { /* check if it is a least squares problem */
    mexErrMsgTxt("Error: size(a,1) < size(a,2).");
  }
  
  x = gsl_matrix_alloc(n, d);
    

  if (has_ab) {
    /* create GSL matrices with the Matlab arrays */
    a = gsl_matrix_alloc(m, n);
    b = gsl_matrix_alloc(m, d);

    m_to_gsl_matrix(a, mxGetPr( prhs[0] ));
    m_to_gsl_matrix(b, mxGetPr( prhs[1] ));
  }


  /* initial approximation prhs[3] */
  if (has_x) {
    /* check dimensions of prhs[3] */
    if ( n != mxGetM(prhs[3]) ) {/* check n */
      mexErrMsgTxt("Error: n ~= size(x,1).");
    }
    if ( d != mxGetN(prhs[3]) ) { /* check d */
      mexErrMsgTxt("Error: d ~= size(x,1).");
    }
    
    if (n+d != n_plus_d) {
      mexErrMsgTxt("Error: n+d ~= total number of columns.");
    }
    
    /* convert x in GSL format and store it in xh */
    m_to_gsl_matrix(x, mxGetPr( prhs[3] ) );
  } 
  
  
  /* Parameter vector p */ 
  if (has_p) {
    if (mxGetN(prhs[5]) != 1 ) { /* check p */
      mexErrMsgTxt("Error: size(p,2) ~= 1.");
    }
    p = NULL;
    
    np = mxGetM(prhs[5]);
//    PRINTF("np = %d", np);
    p = gsl_vector_alloc(np);
    vec_p = gsl_vector_view_array(mxGetPr(prhs[5]), np);
    gsl_vector_memcpy(p, &vec_p.vector);
  } else {
    p = NULL;
  }
 
  

  /*  PRINTF("X\n");
    print_mat(x);*/

  /* optimization options prhs[4] */
  /* default options */
  opt.maxiter = DEF_maxiter;
  opt.epsrel  = DEF_epsrel;
  opt.epsabs  = DEF_epsabs;
  opt.epsgrad = DEF_epsgrad;
  opt.disp    = DEF_disp;
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
	/* lowercase field_name */
	for (c = field_name; *c != '\0'; c++)
	  *c = tolower(*c);
	/* which option */
	IfCheckAndStoreFieldBoundL(maxiter, 1) 
        else IfCheckAndStoreFieldBoundLU(epsrel, 0, 1) 
        else IfCheckAndStoreFieldBoundLU(epsabs, 0, 1) 
        else IfCheckAndStoreFieldBoundLU(epsgrad, 0, 1) 
	else IfCheckAndStoreFieldBoundL(reggamma, 0) 
        else if (! strcmp(field_name, DISP_STR)) {
 	  mxGetString(mxGetFieldByNumber(prhs[4], 0, l), str_buf, STR_MAX_LEN); 
	  /* lowercase str_buf */
	  for (c = str_buf; *c != '\0'; c++) {
	    *c = tolower(*c);
	  }
	    
	  for (i = sizeof(str_disp)/sizeof(str_disp[0]) - 1; i > 0; i--)  {
	    if (!strcmp(str_buf, str_disp[i])) {
	      opt.disp = i;
	      break;
	    }
	  }
	  if (i == 0) {
	    mexWarnMsgTxt("Ignoring optimization option 'disp'. Unrecognized value.");
	  }
 	} else { 
 	  sprintf(err_msg, "Ignoring unrecognized optimization option '%s'.", field_name); 
	  mexWarnMsgTxt(err_msg); 
 	} 
      }
    }
  }
  

  /* --------------- */
  /* Call the solver */
  /* --------------- */
  
  v = gsl_matrix_alloc(n*d, n*d);
  
/*  PRINTF("Hello!! size(a) = (%d,%d), size(b) = (%d,%d), size(x) = (%d,%d), size(v) = (%d,%d), has_p = %d\n", 
      a->size1, a->size2, b->size1, b->size2, x->size1, x->size2, v->size1, v->size2, has_p);*/
  stls(a, b, &s, x, v, &opt, p, has_x, (nlhs > 3) && has_p);
/*  PRINTF("Hello again!!\n");*/

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
/*    temp = mxCreateScalarDouble( opt.fmin ); */
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
  

  if (nlhs > 3 && has_p) {
    plhs[3] = mxCreateDoubleMatrix(np, 1, mxREAL);
    vec_p = gsl_vector_view_array(mxGetPr(plhs[3]), np);
    gsl_vector_memcpy(&vec_p.vector, p);
  }


  /* --------------------- */
  /* Free allocated memory */
  /* --------------------- */

  if (has_p) {
    gsl_vector_free(p);
  }
  
  gsl_matrix_free(a);
  gsl_matrix_free(b);
  gsl_matrix_free(x);
  gsl_matrix_free(v);
}






