#include "slra_mex_utils.h"

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
    MATStoreOption(Mopt, opt, avoid_xi, 0, 1);
  }
}


