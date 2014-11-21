#include "slra.h"

size_t SLRAObject::myObjCnt = 0;
gsl_error_handler_t *SLRAObject::old_gsl_err_h = NULL;

void SLRAObject::myErrorH( const char *reason, const char *F, 
                           int ln, int gsl_err ) {
  throw new Exception("GSL error #%d at %s:%d: %s", ln, 
                      (F != NULL ? F : "<unknown>"),  gsl_err,
                      (reason != NULL ? reason : "<unknown reason>"));
}

SLRAObject::SLRAObject( gsl_vector p_in, gsl_vector ml, gsl_vector nk,
                        gsl_matrix perm, gsl_vector wk, gsl_vector rvec,
                        bool isgcd ) {
  double tmp_n;

  if (old_gsl_err_h != SLRAObject::myErrorH) {
    old_gsl_err_h = gsl_set_error_handler(SLRAObject::myErrorH);
  }

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
    
  myF = new VarproFunction(vecChkNIL(p_in), myS, m-r, matChkNIL(perm), isgcd);
  ++myObjCnt;
}

SLRAObject::~SLRAObject() {
  if (!(--myObjCnt)) {
    Log::deleteLog();
  }
  if (myObjCnt == 0 && old_gsl_err_h != NULL) {
    gsl_set_error_handler(old_gsl_err_h);
    old_gsl_err_h = NULL;
  }
  
  delete myF;
  delete myS;
}

