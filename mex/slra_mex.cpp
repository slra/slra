#include "slra_mex_fun.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  gsl_error_handler_t *new_gsl_err_h = myMexErrorH;
  gsl_error_handler_t *old_gsl_err_h = gsl_set_error_handler(new_gsl_err_h);
  char str_buf[STR_MAX_LEN];
  gsl_matrix rini = { 0, 0, 0, 0, 0, 0 }, rhm = { 0, 0, 0, 0, 0, 0 },
             vhm = { 0, 0, 0, 0, 0, 0 }, psi = { 0, 0, 0, 0, 0, 0 };
  double tmp_n;             
  SLRAObject *slraObj = NULL;
  OptimizationOptions opt;
  
  int was_error = 0;
  try {
    if (nrhs < 3) {     /* Parse required arguments */
      throw new Exception("At least three parameters (p, s, r) are needed.");
    }
    
    slraObj = new SLRAObject(M2vec(prhs[0]), M2vec(mxGetField(prhs[1], 0, ML_STR)), 
                             M2vec(mxGetField(prhs[1], 0, NK_STR)),
                             M2trmat(mxGetField(prhs[1], 0, PERM_STR)),
                             M2vec(mxGetField(prhs[1], 0, WK_STR)), M2vec(prhs[2]) );
    int m =  slraObj->getF()->getNrow();
    int d =  slraObj->getF()->getD();

    /* Parse user supplied options */
    if (nrhs > 3) {
      mexFillOpt(prhs[3], opt, rini, psi, m, m-d); 
    }  
    int mtheta = psi.data != NULL ? psi.size2 : m;

    /* Prepare output info */
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    gsl_vector p_out = M2vec(plhs[0]);
    if (nlhs > 1) {
      mwSize l = 1;
      const char *names[] = { RH_STR, VH_STR, FMIN_STR, ITER_STR, TIME_STR };
      plhs[1] = mxCreateStructArray(1, &l, 5, names);
      mxArray *rh, *vh;
      rhm = M2trmat(rh = mxCreateDoubleMatrix(d, m, mxREAL));
      vhm = M2trmat(vh = mxCreateDoubleMatrix(d*(mtheta-d), d*(mtheta-d), mxREAL));
      mxSetField(plhs[1], 0, RH_STR, rh);
      mxSetField(plhs[1], 0, VH_STR, vh);
    }

    /* Call slra solver and output info */
    slra(slraObj->getF(), &opt, matChkNIL(rini), matChkNIL(psi), 
         vecChkNIL(p_out), matChkNIL(rhm), matChkNIL(vhm));
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
  if (slraObj != NULL) {
    delete slraObj;
  }
  Log::deleteLog();
  if (was_error) {
    mexErrMsgTxt(str_buf);
  }
}
