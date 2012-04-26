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
#include <time.h>

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

typedef struct  {
  gsl_matrix *x;
  gsl_matrix *v;
  gsl_vector *p;
  gsl_matrix *perm;
  opt_and_info opt;
  data_struct s;
  int compdp;
  int given_x;
  slra_opt_data_reshaped params;
  
} slra_r_object; 


SEXP is_slra_object(SEXP ptr) {
  SEXP ans;
  slra_r_object *pSO = NULL;

  PROTECT(ans = allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] = 1;

  /* object is an external pointer */
  if (TYPEOF(ptr) != EXTPTRSXP)
    LOGICAL(ans)[0] = 0;

  /* tag should be 'external matrix' */
  if (LOGICAL(ans)[0] &&
      R_ExternalPtrTag(ptr) != install("SLRA object"))
    LOGICAL(ans)[0] = 0;

  /* pointer itself should not be null */
  if (LOGICAL(ans)[0]) {
    pSO = R_ExternalPtrAddr(ptr);
    if (!pSO)
      LOGICAL(ans)[0] = 0;
  }


  UNPROTECT(1);

  return ans;
}

static void slra_object_finalizer(SEXP ptr) {
  slra_r_object *pSO;

  if (TYPEOF(ptr) != EXTPTRSXP)
    return;

  pSO = R_ExternalPtrAddr(ptr);
  if (!pSO)
    return;
    
  free_memory_reshaped(&pSO->params);

  gsl_vector_free(pSO->p);
  gsl_matrix_free(pSO->x);
  gsl_matrix_free(pSO->v);
  gsl_matrix_free(pSO->perm);
  
  Free(pSO);
  R_ClearExternalPtr(ptr);
}

SEXP create_slra_object(SEXP N, SEXP D, SEXP P, SEXP S, SEXP X, SEXP OPTS, SEXP CDP) {
  slra_r_object *pSO = Calloc(1, slra_r_object);

  SEXP S_A = getListElement(S, "A");
  int *dimS_A = INTEGER(getAttrib(S_A, R_DimSymbol));
  double *s_m = REAL(S_A);
  int *pn = INTEGER(N), *pd = INTEGER(D), *pcompdp = INTEGER(CDP);
  int n = *pn, d = *pd;
  int np;
  gsl_vector_view p_vec;
   
  pSO->compdp =  *pcompdp;
  pSO->given_x = 0;   
  
  /* Variables for input into stls */
  pSO->x = NULL;
  pSO->p = NULL;
  pSO->perm = NULL;
    
  /* Check X */
  pSO->x = gsl_matrix_alloc(n, d);
  if (TYPEOF(X) != NILSXP) {
    pSO->given_x = 1;
    m_to_gsl_matrix(pSO->x, REAL(X));
  }
  pSO->v = gsl_matrix_calloc(n * d, n * d);
  pSO->perm = gsl_matrix_alloc(n + d, n + d);

  /* Copy P vector */
  np = length(P);
  p_vec = gsl_vector_view_array(REAL(P), np);
  pSO->p = gsl_vector_alloc(np);
  gsl_vector_memcpy(pSO->p, &p_vec.vector);

  /* Convert structure specification */
  getScalarListElement(pSO->s.k, S, "k", asInteger, 1);
  slraMatrix2Struct(&pSO->s, s_m, dimS_A[0], dimS_A[1]);

  /* Convert options */
  pSO->opt.disp = getRSLRADispOption(OPTS);
  slraAssignDefOptValue(pSO->opt,method);
  slraAssignDefOptValue(pSO->opt,submethod);
  getRSLRAOption(pSO->opt, OPTS, maxiter, asInteger);
  getRSLRAOption(pSO->opt, OPTS, epsabs, asReal);
  getRSLRAOption(pSO->opt, OPTS, epsrel, asReal);
  getRSLRAOption(pSO->opt, OPTS, epsgrad, asReal);
  getRSLRAOption(pSO->opt, OPTS, epsx, asReal);
  getRSLRAOption(pSO->opt, OPTS, step, asReal);
  getRSLRAOption(pSO->opt, OPTS, tol, asReal);
  getRSLRAOption(pSO->opt, OPTS, reggamma, asReal);
  getRSLRAMethodOption(&pSO->opt, OPTS);
  
  
  /* Make an external pointer envelope */
  SEXP sobj = R_MakeExternalPtr(pSO, install("SLRA object"), R_NilValue);
  R_RegisterCFinalizer(sobj, slra_object_finalizer);

  if (slra_allocate_params(&pSO->params, pSO->p, &pSO->s, pSO->x, pSO->v, 
      &pSO->opt, pSO->given_x, pSO->compdp, pSO->perm, 0) == GSL_EINVAL) {
    error("Error in allocating the params for the computation"); 
  }

  return sobj;
}

slra_r_object *check_slra_object(SEXP ptr) {
  SEXP tchk;

  PROTECT(tchk = is_slra_object(ptr));
  if (LOGICAL(tchk)[0]) {
    return R_ExternalPtrAddr(ptr);
  } else {
    error("pointer provided is not the SLRA object");
  }
  
  return NULL;
}

SEXP optimize_slra(SEXP ptr) {
  slra_r_object *pSO = check_slra_object(ptr);

  if (pSO != NULL) {
    time_t t_b = clock();
    gsl_vector_view x_vec;
    x_vec = gsl_vector_view_array(pSO->x->data, pSO->x->size1 * pSO->x->size2);

    int status = slra_gsl_optimize(&pSO->params, &pSO->opt, &(x_vec.vector), pSO->v);
    pSO->opt.time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;
  }
 

  UNPROTECT(1);
  return NILSXP;
}



SEXP compute_slra_correction(SEXP ptr) {
  SEXP ans = NILSXP;
  slra_r_object *pSO = check_slra_object(ptr);

  if (pSO != NULL) {
    gsl_vector_view x_vec;
    x_vec = gsl_vector_view_array(pSO->x->data, pSO->x->size1 * pSO->x->size2);
    gsl_vector *newp = gsl_vector_alloc(pSO->p->size);  
    gsl_vector_memcpy(newp, pSO->p);
  
    slra_correction_reshaped(newp, &pSO->s, &pSO->params, &(x_vec.vector));

    PROTECT(ans = allocVector(REALSXP, newp->size));
    gsl_vector_view p_vec = gsl_vector_view_array(REAL(ans), newp->size);
    gsl_vector_memcpy(&p_vec.vector, newp);
    UNPROTECT(1);

    gsl_vector_free(newp);
  }

  UNPROTECT(1);
  return ans;
}

SEXP evaluate_slra_f(SEXP ptr, SEXP X) {
  SEXP ans = NILSXP;
  slra_r_object *pSO = check_slra_object(ptr);

  if (pSO != NULL && TYPEOF(X) != NILSXP) {
    int f_len = pSO->params.m * pSO->params.d;
    gsl_vector *x = gsl_vector_alloc(length(X));
    gsl_vector *f = gsl_vector_alloc(f_len);
    gsl_vector_view x_vec = gsl_vector_view_array(REAL(X), length(X));
 
    gsl_vector_memcpy(x, &(x_vec.vector));
    slra_f_reshaped(x, &pSO->params, f);
    
    PROTECT(ans = allocVector(REALSXP, f_len));
    gsl_vector_view f_vec = gsl_vector_view_array(REAL(ans), f_len);
    gsl_vector_memcpy(&f_vec.vector, f);
    UNPROTECT(1);
    
    gsl_vector_free(x);
    gsl_vector_free(f);
  }

  UNPROTECT(1);
  return ans;
}

SEXP evaluate_slra_df(SEXP ptr, SEXP X) {
  SEXP ans = NILSXP;
  slra_r_object *pSO = check_slra_object(ptr);

  if (pSO != NULL && TYPEOF(X) != NILSXP) {
    int f_len = pSO->params.m * pSO->params.d;
    gsl_vector *x = gsl_vector_alloc(length(X));
    gsl_vector *df = gsl_matrix_alloc(f_len,length(X));
    gsl_vector_view x_vec = gsl_vector_view_array(REAL(X), length(X));
 
    gsl_vector_memcpy(x, &(x_vec.vector));
    slra_df_reshaped(x, &pSO->params, df);
    
    PROTECT(ans = allocMatrix(REALSXP, f_len, length(X)));
    gsl_to_m_matrix(REAL(ans), df);
    UNPROTECT(1);
    
    gsl_vector_free(x);
    gsl_matrix_free(df);
  }

  UNPROTECT(1);
  return ans;
}


SEXP return_slra_object(SEXP ptr) {
  SEXP res = NILSXP;
  slra_r_object *pSO = check_slra_object(ptr);

  if (pSO != NULL) {
    SEXP XH, INFO, VXH;
 
    PROTECT(INFO = list3(ScalarInteger(pSO->opt.iter), ScalarReal(pSO->opt.time), 
                         ScalarReal(pSO->opt.fmin)));
    SET_TAG(INFO, install("iter"));
    SET_TAG(CDR(INFO), install("time"));
    SET_TAG(CDDR(INFO), install("fmin"));
    PROTECT(XH = allocMatrix(REALSXP, pSO->x->size1, pSO->x->size2));
    PROTECT(VXH = allocMatrix(REALSXP, pSO->v->size1, pSO->v->size2));
    gsl_to_m_matrix(REAL(XH), pSO->x); 
    gsl_to_m_matrix(REAL(VXH), pSO->v); 

    PROTECT(res = list3(XH, INFO, VXH));
    SET_TAG(res, install("xh"));
    SET_TAG(CDR(res), install("info"));
    SET_TAG(CDDR(res), install("vxh"));

    UNPROTECT(4);
  }

  UNPROTECT(1);
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



