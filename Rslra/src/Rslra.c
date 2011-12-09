#include <stdarg.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <R_ext/Print.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/BLAS.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>

#include "slra/slra.h"

/* default constants for the exit condition */
static SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = GET_NAMES(list);
  int i;

  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

#define getScalarListElement(trg, list, str, coerce, def)         \
  do {                                                            \
    SEXP __tmp = getListElement(list, str);                       \
    trg = (__tmp != R_NilValue ? coerce(__tmp) : (def));          \
  } while(0)

#define getRSLRAOption(opt, list, field, coerce)          \
  do {                                                    \
    getScalarListElement(opt.field, list, #field, coerce, \
        SLRA_DEF_##field);  \
  } while(0)

static int getRSLRADispOption( SEXP OPTS ) {
  SEXP str_value_sexp;
  
  if (TYPEOF((str_value_sexp = getListElement(OPTS, "disp"))) == STRSXP) {
    return slraString2Disp(CHAR(STRING_ELT(str_value_sexp, 0)));    
  }
  
  return SLRA_DEF_disp;
}

static void getRSLRAMethodOption( opt_and_info *popt, SEXP OPTS ) {
  SEXP str_val_sexp;

  if (TYPEOF((str_val_sexp = getListElement(OPTS, "method"))) == STRSXP) {
    slraString2Method(CHAR(STRING_ELT(str_val_sexp, 0)), popt);
  } else {
    slraAssignDefOptValue((*popt), method);
    slraAssignDefOptValue((*popt), submethod);
  }
}

/********************************
 * C wrapper function.
 * Assumed that A and B are matrices and all dimensions are checked.
 * It is assumed that X exists.
 * Assumes that S is a list
 ******************************/
SEXP rslra(SEXP N, SEXP D, SEXP P, SEXP S, SEXP X, SEXP OPTS, SEXP CDP) {
  static const char *str_codes = " THUE";
  SEXP S_A = getListElement(S, "A");
  int *dimS_A = INTEGER(getAttrib(S_A, R_DimSymbol));
  double *s_m = REAL(S_A);
  int *pn = INTEGER(N), *pd = INTEGER(D), *pcompdp = INTEGER(CDP);
  int n = *pn, d = *pd, compdp = *pcompdp;
  int np;
  gsl_vector_view p_vec;
   
  /* Variables for input into stls */
  gsl_matrix *x = NULL, *v;
  gsl_vector *p = NULL;
  gsl_matrix *perm = NULL;
  data_struct s;
  opt_and_info opt;

  /* Check X */
  x = gsl_matrix_alloc(n, d);
  if (TYPEOF(X) != NILSXP) {
    m_to_gsl_matrix(x, REAL(X));
  }
  v = gsl_matrix_calloc(n * d, n * d);
  perm = gsl_matrix_alloc(n + d, n + d);

  /* Copy P vector */
  np = length(P);
  p_vec = gsl_vector_view_array(REAL(P), np);
  p = gsl_vector_alloc(np);
  gsl_vector_memcpy(p, &p_vec.vector);

  /* Convert structure specification */
  getScalarListElement(s.k, S, "k", asInteger, 1);
  slraMatrix2Struct(&s, s_m, dimS_A[0], dimS_A[1]);

  /* Convert options */
  opt.disp = getRSLRADispOption(OPTS);
  slraAssignDefOptValue(opt,method);
  slraAssignDefOptValue(opt,submethod);
  getRSLRAOption(opt, OPTS, maxiter, asInteger);
  getRSLRAOption(opt, OPTS, epsabs, asReal);
  getRSLRAOption(opt, OPTS, epsrel, asReal);
  getRSLRAOption(opt, OPTS, epsgrad, asReal);
  getRSLRAOption(opt, OPTS, epsx, asReal);
  getRSLRAOption(opt, OPTS, step, asReal);
  getRSLRAOption(opt, OPTS, tol, asReal);
  getRSLRAOption(opt, OPTS, reggamma, asReal);
 
  slra(p, &s, x, v, &opt, TYPEOF(X) != NILSXP, compdp, perm, 0); 

  /* Form the result */
  SEXP XH, INFO, VXH, PH, res;

  PROTECT(INFO = list3(ScalarInteger(opt.iter), ScalarReal(opt.time), 
                       ScalarReal(opt.fmin)));
  SET_TAG(INFO, install("iter"));
  SET_TAG(CDR(INFO), install("time"));
  SET_TAG(CDDR(INFO), install("fmin"));
  PROTECT(XH = allocMatrix(REALSXP, n, d));
  PROTECT(VXH = allocMatrix(REALSXP, n*d, n*d));
  gsl_to_m_matrix(REAL(XH), x); 
  gsl_to_m_matrix(REAL(VXH), v); 

  if (!compdp) {
    PROTECT(res = list3(XH, INFO, VXH));
    SET_TAG(res, install("xh"));
    SET_TAG(CDR(res), install("info"));
    SET_TAG(CDDR(res), install("vxh"));

    UNPROTECT(4);
  } else {
    PROTECT(PH = allocVector(REALSXP, np));

    p_vec = gsl_vector_view_array(REAL(PH), np);
    gsl_vector_memcpy(&p_vec.vector, p);
  
    PROTECT(res = list4(XH, INFO, VXH, PH));
    SET_TAG(res, install("xh"));
    SET_TAG(CDR(res), install("info"));
    SET_TAG(CDDR(res), install("vxh"));
    SET_TAG(CDR(CDDR(res)), install("ph"));

    UNPROTECT(5);
  }

  gsl_vector_free(p);
  gsl_matrix_free(x);
  gsl_matrix_free(v);
  gsl_matrix_free(perm);

  return res;
}


SEXP rtls(SEXP A, SEXP B) {
  int *dimA = INTEGER(getAttrib(A, R_DimSymbol));
  int *dimB = INTEGER(getAttrib(B, R_DimSymbol));
  SEXP X;

  /* Variables for input into stls */
  gsl_matrix *a, *b, *x;
  
  m_to_gsl_matrix(a = gsl_matrix_alloc(dimA[0], dimA[1]), REAL(A));
  m_to_gsl_matrix(b = gsl_matrix_alloc(dimB[0], dimB[1]), REAL(B));
  x = gsl_matrix_calloc(dimA[1], dimB[1]);

  tls(a, b, x);

  PROTECT(X = allocMatrix(REALSXP, dimA[1], dimB[1]));

  gsl_to_m_matrix(REAL(X), x); 

  UNPROTECT(1);

  gsl_matrix_free(a);
  gsl_matrix_free(b);
  gsl_matrix_free(x);

  return X;
}



