#include "slra_mex_fun.h"
#include <typeinfo>
using namespace std;
#include "class_handle.hpp"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  char str_buf[STR_MAX_LEN];
  int was_error = 0;
  try {
    if (nrhs < 1 || mxGetString(prhs[0], str_buf, STR_MAX_LEN)) {
      throw new Exception("The first argument should be a command string.");
    }
    
    if (!strcmp("new", str_buf)) { /* Create new SLRA object */
      if (nrhs < 4) {     
        throw new Exception("At least three parameters (p, s, r) are needed.");
      }
      const mxArray *as = prhs[2];
      gsl_vector isgcd = M2vec(mxGetField(as, 0, GCD_STR));
      SLRAObject *slraObj = new SLRAObject(M2vec(prhs[1]), 
         M2vec(mxGetField(as, 0, ML_STR)), M2vec(mxGetField(as, 0, NK_STR)),
         M2trmat(mxGetField(as, 0, PERM_STR)), M2vec(mxGetField(as, 0, WK_STR)), 
         M2vec(prhs[3]), (isgcd.data != NULL) && (*isgcd.data));
                               
      plhs[0] = convertPtr2Mat<SLRAObject>(slraObj);                             
      return;
    }
    
    if (nrhs < 2) {
      throw new Exception("The second argument should be an SLRA object.");
    }
    
    if (!strcmp("delete", str_buf)) { /* Delete an SLRA object */
      destroyObject<SLRAObject>(prhs[1]);
      return;
    }

    SLRAObject * slraObj = convertMat2Ptr<SLRAObject>(prhs[1]);
    int m =  slraObj->getF()->getNrow(), d =  slraObj->getF()->getD();
    

    if (!strcmp("optimize", str_buf)) {
      OptimizationOptions opt;
      gsl_matrix rini = { 0, 0, 0, 0, 0, 0 }, psi = { 0, 0, 0, 0, 0, 0 }, 
                 rhm = { 0, 0, 0, 0, 0, 0 }, vhm =  { 0, 0, 0, 0, 0, 0 };
      if (nrhs > 2) {
        mexFillOpt(prhs[2], opt, rini, psi, m, m-d); 
      }  
      int mtheta = psi.data != NULL ? psi.size2 : m;
      /* Prepare output info */
      plhs[0] = mxCreateDoubleMatrix(slraObj->getF()->getNp(), 1, mxREAL);
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
      return;
    }

    if (nlhs <= 0) {
      throw new Exception("Output arguments should be provided.");        
    }

    if (!strcmp("getRini", str_buf)) {
      gsl_matrix rini = M2trmat(plhs[0] = mxCreateDoubleMatrix(d, m, mxREAL));
      slraObj->getF()->computeDefaultRTheta(&rini);
      return;
    }
    if (!strcmp("getM", str_buf)) {
      plhs[0] = mxCreateDoubleScalar(slraObj->getF()->getNrow());
      return;
    }

    if (nrhs < 3) {
      throw new Exception("The third argument should be R.");        
    }
    gsl_matrix R = M2trmat(prhs[2]);
    
    if (!strcmp("func", str_buf)) {
      double res;
      slraObj->getF()->computeFuncAndGrad(&R, &res, NULL, NULL);      
      plhs[0] = mxCreateDoubleScalar(res);
      return;
    }
    if (!strcmp("grad", str_buf)) {
      gsl_matrix grad_m = M2trmat(plhs[0] = mxCreateDoubleMatrix(d, m, mxREAL));
      slraObj->getF()->computeFuncAndGrad(&R, NULL, NULL, &grad_m);      

      return;
    }
    if (!strcmp("getPh", str_buf)) {
      gsl_vector pout = 
           M2vec(plhs[0] = mxCreateDoubleMatrix(slraObj->getF()->getNp(), 1, mxREAL));
      slraObj->getF()->computePhat(&pout, &R);      
      return;
    }
  } catch (Exception *e) {
    strncpy0(str_buf, e->getMessage(), STR_MAX_LEN);
    was_error = 1;
    delete e;
  } 

  if (was_error) {
#ifdef   BUILD_MEX_WINDOWS
    mexPrintf("Error using slra MEX object:\n");
    mexPrintf("%s\n", str_buf);
#else /* BUILD_MEX_WINDOWS */
    mexErrMsgTxt(str_buf);
#endif /* BUILD_MEX_WINDOWS */      
  }
}
