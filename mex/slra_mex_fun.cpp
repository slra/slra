#include "slra_mex_fun.h"

void myMexErrorH( const char *reason, const char *F, int ln, int gsl_err ) {
  throw new Exception("GSL error #%d at %s:%d: %s", ln, 
                      (F != NULL ? F : "<unknown>"),  gsl_err,
                      (reason != NULL ? reason : "<unknown reason>"));
}

gsl_matrix M2trmat( const mxArray * mat ) {
  gsl_matrix res = { 0, 0, 0, 0, 0, 0 };
  if (mat != NULL && mxGetN(mat) != 0 && mxGetM(mat) != 0) {
    res = gsl_matrix_const_view_array(mxGetPr(mat),mxGetN(mat),mxGetM(mat)).matrix;
  }
  return res;
}

gsl_vector M2vec( const mxArray * mat ) {
  gsl_vector res = { 0, 0, 0, 0, 0 };
  if (mat != NULL && mxGetN(mat) != 0 && mxGetM(mat) != 0) {
    res = gsl_vector_const_view_array(mxGetPr(mat), 
                                      mxGetN(mat) * mxGetM(mat)).vector;
  }
  return res;
}

char *M2Str( mxArray *myMat, char *str, size_t max_len ) {
  if (myMat == NULL) {
    *str = 0;
  } else {   
    mxGetString(myMat, str, max_len); 
  }  
  return str;
}

void mexFillOpt( const mxArray *Mopt, OptimizationOptions &opt, 
                 gsl_matrix & Rini, gsl_matrix &Psi, size_t m, size_t r ) {
  char str_buf[STR_MAX_LEN];
  if (! mxIsStruct(Mopt)) {
    mexWarnMsgTxt("Ignoring 'opt'. The optimization options "
	          "should be passed in a structure.");
  } else {
    Log::str2DispLevel(M2Str(mxGetField(Mopt, 0, DISP_STR), str_buf, 
                             STR_MAX_LEN));
    Rini = M2trmat(mxGetField(Mopt, 0, RINI_STR));
    if (Rini.data != NULL && (Rini.size2 != (m-r) || Rini.size1 != m)) {
      throw new Exception("Incorrect Rini\n");   
    }

    Psi = M2trmat(mxGetField(Mopt, 0, PSI_STR));
    if (Psi.data != NULL) {
      if (Psi.size1 != m || Psi.size2 == 0 || Psi.size2 > m) {
        throw new Exception("Incorrect Psi\n");   
      }
      
      if ((m-r) >= Psi.size2) {
        throw new Exception("Rank reduction and psi incompatible\n");   
      }
    } 
    opt.str2Method(M2Str(mxGetField(Mopt, 0, METHOD_STR), str_buf, 
                                   STR_MAX_LEN));
    MATStoreOption(Mopt, opt, maxiter, 0, numeric_limits<int>::max());
    MATStoreOption(Mopt, opt, epsabs, 0, 1);
    MATStoreOption(Mopt, opt, epsrel, 0, 1);
    MATStoreOption(Mopt, opt, epsgrad, 0, 1);
    MATStoreOption(Mopt, opt, epsx, 0, 1);
    MATStoreOption(Mopt, opt, maxx, 0, numeric_limits<double>::max());
    MATStoreOption(Mopt, opt, step, 0, 1);
    MATStoreOption(Mopt, opt, tol, 0, 1);
    MATStoreOption(Mopt, opt, reggamma, 0, numeric_limits<double>::max());
    MATStoreOption(Mopt, opt, ls_correction, 0, 1);
  }
}

size_t SLRAObject::myObjCnt = 0;


SLRAObject::SLRAObject( gsl_vector p_in, gsl_vector ml, gsl_vector nk,
                        gsl_matrix perm, gsl_vector wk, gsl_vector rvec ) :
                        old_gsl_err_h(NULL)  {
  double tmp_n;
  if (ml.size == 0) {
    throw new Exception("s.m should be a nonempty vector");   
  }
    
  if (nk.size == 0) {
    tmp_n = compute_n(&ml, p_in.size);
    nk = gsl_vector_view_array(&tmp_n, 1).vector;
  }    
  size_t np = compute_np(&ml, &nk);

  if (p_in.size < np) {
    throw new Exception("Size of vector p less than needed");   
  } else if (p_in.size > np) {
    throw new Exception("Size of vector p exceeds structure requirements");   
  } 

  myS = createMosaicStructure(&ml, &nk, vecChkNIL(wk));
  size_t m = (perm.data == NULL ? myS->getM() : perm.size2);
  int r = (rvec.size == 0 ? m - 1 : gsl_vector_get(&rvec, 0));
  
  if (r <= 0 || r >= m) {
    throw new Exception("Incorrect rank\n");   
  }
    
  myF = new CostFunction(vecChkNIL(p_in), myS, m-r, matChkNIL(perm));
  old_gsl_err_h = gsl_set_error_handler(myMexErrorH);
  ++myObjCnt;
}

SLRAObject::~SLRAObject() {
  if (!(--myObjCnt)) {
    Log::deleteLog();
  }
  if (old_gsl_err_h != myMexErrorH) {
    gsl_set_error_handler(old_gsl_err_h);
  }
  
  delete myF;
  delete myS;
}

