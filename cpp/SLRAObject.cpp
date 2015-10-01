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

  myS = createMosaicStructure(&ml, &nk, vecChkNIL(wk), matChkNIL(perm));
  size_t m = (perm.data == NULL ? myS->getM() : perm.size2);
  int r = (rvec.size == 0 ? m - 1 : gsl_vector_get(&rvec, 0));
  
  if (r <= 0 || r >= m) {
    throw new Exception("Incorrect rank\n");   
  }
    
  myF = new VarproFunction(vecChkNIL(p_in), myS, m-r, NULL, isgcd);
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

void SLRAObject::optimize( OptimizationOptions* opt, gsl_matrix *Rini,
           gsl_matrix *Psi, gsl_vector *p_out, gsl_matrix *r_out, 
           gsl_matrix *v_out, gsl_matrix *Rs, gsl_matrix *info ) { 
  NLSVarpro *optFun = NULL;
  gsl_vector *x = NULL;
  double old_reg = myF->getReggamma();
  
  try { 
    time_t t_b = clock();

    myF->setReggamma(opt->reggamma);
    if (Psi != NULL && Psi->size1 != myF->getNrow()) {
      opt->avoid_xi = 1;
    }
    
    if (myF->isGCD()) {
      opt->ls_correction = 1;
    }
    if (opt->ls_correction) {
      if (opt->avoid_xi) {
        if (Psi != NULL) {
          optFun = new NLSVarproPsiVecRCorrection(*myF, Psi);
        } else {
          optFun = new NLSVarproVecRCorrection(*myF);
        }
      } else {
        optFun = new NLSVarproPsiXICorrection(*myF, Psi);
      }
    } else {
      if (opt->avoid_xi) {
        if (Psi != NULL) {
          optFun = new NLSVarproPsiVecRCholesky(*myF, Psi);
        } else {
          optFun = new NLSVarproVecRCholesky(*myF);
        }
      } else {
        optFun = new NLSVarproPsiXICholesky(*myF, Psi);
      }
    }
    x = gsl_vector_alloc(optFun->getNvar());

    MyIterationLogger itLog(optFun, Rs, info);

    if (Rini == NULL) {  
      Log::lprintf(Log::LOG_LEVEL_ITER, 
           "R not given - computing initial approximation.\n");    
      optFun->computeDefaultx(x);
    } else {
      optFun->RTheta2x(Rini, x);
    }
    if (opt->avoid_xi) {
      opt->method = SLRA_OPT_METHOD_LMPINV;
    }

    if (opt->method == SLRA_OPT_METHOD_LMPINV) {
      opt->lmpinvOptimize(optFun, x, &itLog);
    } else {
      opt->gslOptimize(optFun, x, v_out, &itLog);
    } 

    if (p_out != NULL) {
      optFun->computePhat(p_out, x);
    }
    opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;
    if (r_out != NULL) {
      optFun->x2RTheta(r_out, x);
    }

    throw (Exception *)NULL; /* Throw NULL exception to unify deallocation */
  } catch ( Exception *e ) {
    if (optFun != NULL)  {
      delete optFun;
    }
    gsl_vector_free_ifnull(x);
    
    if (e != NULL) { /* Abnormal termination only if e is normal exception */
      throw;  
    }
    myF->setReggamma(old_reg);
  }
}