
#if defined(BUILD_R_PACKAGE)

#include <R.h>
#define PRINTF Rprintf
#define MYWARNING Rprintf

#elif defined(BUILD_MEX_OCTAVE) || defined(BUILD_MEX_MATLAB)

#include "mex.h"
#define PRINTF mexPrintf
#define MYWARNING mexWarnMsgTxt

#else

#include <stdio.h>
#define PRINTF printf
#define MYWARNING printf

#endif

#define DEBUGINT(x) PRINTF("%s = %d\n", #x, x)
#define DEBUGDOUBLE(x) PRINTF("%s = %f\n", #x, x)

/* field names for s */
#define ML_STR "m"
#define NK_STR "n"
#define PERM_STR "phi"
#define WK_STR "w"

/* field names for opt */
#define RINI_STR "Rini"
#define DISP_STR "disp"
#define METHOD_STR "method"

/* names for output */
#define PH_STR "ph"
#define INFO_STR "info"


/* field names for info */
#define RH_STR "Rh"
#define VH_STR "Vh"
#define FMIN_STR "fmin"
#define ITER_STR "iter"
#define TIME_STR "time"
#define PSI_STR "psi"
          
#define mymax(a, b) ((a) > (b) ? (a) : (b)) 
#define mymin(a, b) ((a) < (b) ? (a) : (b))
         
/*
 * tmv_prod_new: block-Toeplitz banded matrix p =  T * v
 * T - storage for [t_s-1' ... t_1' t_0 t_1 ... t_s-1].
 * m = number of block rows / columns
 */ 
void tmv_prod_vector( gsl_vector *T, int s, gsl_vector* v, int m, 
         gsl_vector* p );
void tmv_prod_new( gsl_matrix *T, int s,  gsl_vector *v, int m, gsl_vector *p, 
         double beta = 0.0 );

void ls_solve( const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *X ); 
         
void copyLowerTrg( gsl_matrix * dest, const gsl_matrix *src  );         
void shiftLowerTrg( gsl_matrix * dest, const gsl_matrix *src  );         

void print_mat(const gsl_matrix*);
void print_mat_tr(const gsl_matrix*);
void print_arr(const double*, int);

void print_vec(const gsl_vector*);

int read_mat( gsl_matrix *a, const char * filename );
int read_vec( gsl_vector *a, const char * filename );
int read_vec_uint( gsl_vector_uint *a, const char * filename );

void gsl_matrix_vectorize(double*, const gsl_matrix*);
void gsl_matrix_vec_inv(gsl_matrix*, const double*);


               
int compute_np( gsl_vector* ml, gsl_vector *nk );               
int compute_n( gsl_vector* ml, int np );               


const gsl_vector *vecChkNIL( const gsl_vector &vec );
gsl_vector *vecChkNIL( gsl_vector &vec  );
gsl_matrix *matChkNIL( gsl_matrix &mat_vw );

void tolowerstr( char * str );


/** Class for logging */
class Log {
public:  
  enum Level {
    LOG_LEVEL_OFF = 0,
    LOG_LEVEL_FINAL,
    LOG_LEVEL_NOTIFY,
    LOG_LEVEL_ITER
  };

  static void lprintf( char *format, ... );  
  static void lprintf( Level level, char *format, ... );  

  /** Initialize disp field from string */
  static void str2DispLevel( const char *str );
  static void setMaxLevel( Level maxLevel );
  static Level getMaxLevel();

  static void deleteLog();
  
private:
  static const int MSG_MAX = 200;
  char myMsg[MSG_MAX];
  Level myMaxLevel;

  static Log *getLog();


  Log() : myMaxLevel(LOG_LEVEL_NOTIFY) {};
  ~Log() {};
  static Log *myLogInstance;
};

