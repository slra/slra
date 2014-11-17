/* stls.c: implementations of the functions from stls.h */
#include <time.h>
#include "slra.h"
#include "Timer.h"

class MyIterationLogger : public IterationLogger {
  Timer myTimer;
  double myStartTime;
  NLSVarpro *myFun;
  gsl_matrix *myRs;
  gsl_matrix *R;
  gsl_matrix *myInfo;
public:
  MyIterationLogger( NLSVarpro *fun, gsl_matrix *Rs, gsl_matrix *info ) {
    myFun = fun;
    myRs = Rs;
    myInfo = info;  
    myTimer.start();    
  }
  virtual void reportIteration( int no, const gsl_vector *x, double fmin, 
                                   const gsl_vector *grad ) {
    double tm = myTimer.getElapsedTime();
    double grad_nrm = (grad != NULL) ? grad_nrm = gsl_blas_dnrm2(grad) : NULL;
                                   
    if (Log::getMaxLevel() >= Log::LOG_LEVEL_ITER) {
      double x_norm = gsl_blas_dnrm2(x);
      Log::lprintf("%3u: f0 = %15.10e, ||x|| = %10.8f", no, fmin, x_norm);
      if (grad != NULL) {
        Log::lprintf(", ||f0'|| = %15.7e", grad_nrm);
      }
      Log::lprintf("\n");
    }  
    if (myRs != NULL && no >= 0 && no < myRs->size1) {
      gsl_vector RsRow = gsl_matrix_row(myRs, no).vector;
      gsl_matrix R = gsl_matrix_view_vector(&RsRow, myFun->getM(), 
                                            myFun->getD()).matrix;
      myFun->x2RTheta(&R, x);                                            
    }
    
    if (myInfo != NULL && no >= 0 && no < myInfo->size1) {
      gsl_vector InfoRow = gsl_matrix_row(myInfo, no).vector;
      gsl_vector_set(&InfoRow, 0, tm);
      if (InfoRow.size > 1) {
        gsl_vector_set(&InfoRow, 1, fmin);
      }
      if (InfoRow.size > 2 && grad != NULL) {
        gsl_vector_set(&InfoRow, 2, grad_nrm);
      }
    }
  }
};

void slra( VarproFunction *costFun, OptimizationOptions* opt, gsl_matrix *Rini, 
           gsl_matrix *Psi, gsl_vector *p_out, gsl_matrix *r_out, 
           gsl_matrix *v_out, gsl_matrix *Rs, gsl_matrix *info ) { 
  NLSVarpro *optFun = NULL;
  gsl_vector *x = NULL;
  double old_reg = costFun->getReggamma();
  
  try { 
    time_t t_b = clock();

    costFun->setReggamma(opt->reggamma);
    if (costFun->isGCD()) {
      opt->ls_correction = 1;
    }
    if (opt->ls_correction) {
      if (opt->avoid_xi) {
        optFun = new NLSVarproVecRCorrection(*costFun);
      } else {
        optFun = new NLSVarproPsiXICorrection(*costFun, Psi);
      }
    } else {
      if (opt->avoid_xi) {
        optFun = new NLSVarproVecRCholesky(*costFun);
      } else {
        optFun = new NLSVarproPsiXICholesky(*costFun, Psi);
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
    costFun->setReggamma(old_reg);
  }
}
