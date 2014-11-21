#include <time.h>
#include "Timer.h"

class MyIterationLogger : public IterationLogger {
  Timer myTimer;
  double myStartTime;
  NLSVarpro *myFun;
  gsl_matrix *myRs;
  gsl_matrix *R;
  gsl_matrix *myInfo;
public:
  MyIterationLogger( NLSVarpro *fun, gsl_matrix *Rs, gsl_matrix *info ) :
    myFun(fun), myRs(Rs), myInfo(info) {
    myTimer.start();
  }
  virtual void reportIteration( int no, const gsl_vector *x, double fmin, 
                               const gsl_vector *grad );
};


