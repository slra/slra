#include "slra.h"

class SLRAObject {
  Structure *myS;
  VarproFunction *myF;
  static void myErrorH( const char *reason, const char *F, int ln, int gsl_err );
  static gsl_error_handler_t *old_gsl_err_h;
  static size_t myObjCnt;
public:
  SLRAObject( gsl_vector p_in, gsl_vector ml, gsl_vector nk,
              gsl_matrix perm, gsl_vector wk, gsl_vector rvec,
              bool isgcd = false );
  virtual ~SLRAObject();
    
  Structure *getS() { return myS; }
  VarproFunction *getF() { return myF; }
};

