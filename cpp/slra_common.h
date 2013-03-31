
#if defined(BUILD_R_PACKAGE)

#include <R.h>
#define PRINTF Rprintf
#define MYWARNING Rprintf
#define FLUSH()

#elif defined(BUILD_MEX_OCTAVE) || defined(MATLAB_MEX_FILE)

#include "mex.h"
#define PRINTF mexPrintf
#define FLUSH() mexEvalString("drawnow;")
#define MYWARNING mexWarnMsgTxt

#else

#include <stdio.h>
#define PRINTF printf
#define FLUSH()
#define MYWARNING printf

#endif

#define DEBUGINT(x) PRINTF("%s = %d\n", #x, x)
#define DEBUGDOUBLE(x) PRINTF("%s = %f\n", #x, x)

/* field names for s */
#define ML_STR "m"
#define NK_STR "n"
#define PERM_STR "phi"
#define WK_STR "w"
#define GCD_STR "gcd"

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
         
#define gsl_matrix_free_ifnull(M)    if (M != NULL) gsl_matrix_free(M)
#define gsl_vector_free_ifnull(V)    if (V != NULL) gsl_vector_free(V)
     
     
/** Function that creates appropriate Mosaic structure from weights and checks its
 *    consistence with \f$n_p\f$
 * @ingroup MainFunctions 
 * @param [in]     ml      Vector of \f$\bf m\f$ \sa MosaicHStructure
 * @param [in]     nk      Vector of \f$\bf n\f$ \sa MosaicHStructure  
 * @param [in]     wk      Vector of weights. If NULL, or of sizes \f$q\f$, \f$qN\f$ - 
 *                         MosaicHStructure is constructed, otherwise - WMosaicHStructure
 * @param [in]     d       Rank reduction  (\f$m-r\f$) 
 */                
Structure *createMosaicStructure( gsl_vector * ml,  gsl_vector *nk, 
               gsl_vector * wk );

         
/*
 * tmv_prod_new: block-Toeplitz banded matrix p =  T * v
 * T - storage for [t_s-1' ... t_1' t_0 t_1 ... t_s-1].
 * m = number of block rows / columns
 */ 
void tmv_prod_vector( gsl_vector *T, size_t s, const gsl_vector* v, size_t m, 
         gsl_vector* p );
void tmv_prod_new( gsl_matrix *T, size_t s, const gsl_vector *v, size_t m, gsl_vector *p, 
         double beta = 0.0 );

void ls_solve( const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *X ); 
         
void copyLowerTrg( gsl_matrix * dest, const gsl_matrix *src  );         
void shiftLowerTrg( gsl_matrix * dest, const gsl_matrix *src  );         

void print_mat(const gsl_matrix*);
void print_mat_tr(const gsl_matrix*);
void print_arr(const double*, size_t);

void print_vec(const gsl_vector*);

int read_mat( gsl_matrix *a, const char * filename );
int read_vec( gsl_vector *a, const char * filename );
int read_vec_uint( gsl_vector_uint *a, const char * filename );

void gsl_matrix_vectorize(double*, const gsl_matrix*);
void gsl_matrix_vec_inv(gsl_matrix*, const double*);


               
size_t compute_np( gsl_vector* ml, gsl_vector *nk );               
size_t compute_n( gsl_vector* ml, size_t np );               


const gsl_vector *vecChkNIL( const gsl_vector &vec );
gsl_vector *vecChkNIL( gsl_vector &vec  );
gsl_matrix *matChkNIL( gsl_matrix &mat_vw );

void tolowerstr( char * str );

char *strncpy0( char *dest, const char * src, size_t buf_len );


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
  static const size_t MSG_MAX = 200;
  char myMsg[MSG_MAX];
  Level myMaxLevel;

  static Log *getLog();


  Log() : myMaxLevel(LOG_LEVEL_NOTIFY) {};
  ~Log() {};
  static Log *myLogInstance;
};

