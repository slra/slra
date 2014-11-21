#include "slra.h"

void MyIterationLogger ::reportIteration( int no, const gsl_vector *x, double fmin, 
                                   const gsl_vector *grad ) {
    double tm = myTimer.getElapsedTime();
    double grad_nrm = (grad != NULL) ? gsl_blas_dnrm2(grad) : 0;
                                   
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



