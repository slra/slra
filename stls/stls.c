/* stls.c: implementations of the functions from stls.h */

#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include "stls.h"


#define mymax(a,b) ((a) > (b) ? (a) : (b)) 
#define mymin(a,b) ((a) < (b) ? (a) : (b))

/*********************/
/* stls: STLS solver */
/*********************/
/* on exit: zero success, otherwise error code: EITER, GSL_ETOLF, GSL_ETOLX, GSL_ETOLG */

int stls(gsl_matrix* a, gsl_matrix* b, const data_struct* s, 
         gsl_matrix* x, gsl_matrix* v, opt_and_info* opt)
{
  int m, n, d, status;
  char method = 'l';
  int status_dx, status_grad, k;
  double g_norm;
  w_data w;
  stls_opt_data params;
  gsl_vector_view x_vec;
  time_t t_b;

  /* LM */
  const gsl_multifit_fdfsolver_type *Tlm
    = gsl_multifit_fdfsolver_lmder;
  gsl_multifit_fdfsolver* solverlm;
  gsl_multifit_function_fdf fdflm;
  gsl_vector *g;

  /* QN */
  double stepqn = 0.001; /* ??? */
  const gsl_multimin_fdfminimizer_type *Tqn
    = gsl_multimin_fdfminimizer_vector_bfgs;
  gsl_multimin_fdfminimizer* solverqn;
  gsl_multimin_function_fdf fdfqn;

  /* NM */
  double size;
  gsl_vector *stepnm;
  const gsl_multimin_fminimizer_type *Tnm
    = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer* solvernm;
  gsl_multimin_function fnm;


  t_b = clock();

  /* constants */
  m = a->size1;
  n = a->size2;
  d = b->size2;

  /* find Wk */
  s2w(s, &w);

  /* set the parameters */
  params.m = m;
  params.n = n;
  params.d = d;


  params.a = a;
  params.b = b;
  params.w = &w;
  params.k = s->k;  
  params.n_plus_d = n + d;   
  params.n_times_d = n * d;   
  params.k_times_d = s->k * d; 
  params.k_times_d_times_s = params.k_times_d * w.s;
  params.k_times_d_times_s_minus_1 = params.k_times_d_times_s - 1;
  params.m_times_d = m * d;
  params.m_div_k = (int) m / s->k;
  params.s_minus_1 = w.s - 1;



  params.one = 1;

  
  /* Preallocate memory for f and df */
  params.x_ext = gsl_matrix_calloc(params.w->a[0]->size1, params.k_times_d);
  params.yr = gsl_vector_alloc(params.m_times_d);  

  /*CholGam */
  params.rb = (double*) calloc(params.m_times_d * params.k_times_d_times_s, sizeof(double));
  params.ldwork = 1 + (params.s_minus_1 + 1)* params.k_times_d * params.k_times_d +  /* pDW */ 
                       3 * params.k_times_d + /* 3 * K */
                       mymax(params.s_minus_1 + 1,  params.m_div_k - 1 - params.s_minus_1) * params.k_times_d * params.k_times_d; /* Space needed for MB02CV */

  params.dwork  = (double*) malloc((size_t)params.ldwork * sizeof(double));
  params.gamma = gsl_matrix_alloc(params.k_times_d, params.k_times_d_times_s);  
  params.gamma_vec = (double*) malloc(params.k_times_d * params.k_times_d_times_s *sizeof(double));

  params.tmp   = gsl_matrix_alloc(params.k_times_d, params.w->a[0]->size1);



  /* New CholGam */
  params.bx_ext =  gsl_matrix_alloc(params.n_plus_d, d);

  params.rb2 = (double*) malloc(params.m_times_d * params.k_times_d_times_s * sizeof(double));
  params.brg_rb = (double*) malloc(params.m_div_k * d * params.d_times_s * sizeof(double));
  params.d_times_s = d * w.s;
  params.d_times_m_div_k = d* (int) m / s->k;
  params.d_times_s_minus_1 = params.d_times_s - 1;

  params.brg_ldwork = 1 + (params.s_minus_1 + 1)* d * d +  /* pDW */ 
                       3 * d + /* 3 * K */
                       mymax(params.s_minus_1 + 1,  params.m_div_k - 1 - params.s_minus_1) * d * d; /* Space needed for MB02CV */

  params.brg_gamma_vec = (double*) malloc(d * params.d_times_s * sizeof(double));
  params.brg_gamma = gsl_matrix_alloc(d, params.d_times_s);
  params.brg_dwork  = (double*) malloc((size_t)params.brg_ldwork * sizeof(double));
  params.brg_tmp   = gsl_matrix_alloc(d, params.n_plus_d);

  params.brg_yr = gsl_vector_alloc(params.m_times_d);  
  params.brg_f = gsl_vector_alloc(params.m_times_d);  

  
  
  /* Jacobian*/
  params.jres1 = malloc(params.m_times_d * params.n_times_d * sizeof(double));
  params.jres2  = malloc( params.m_times_d * sizeof(double));

  params.dgamma = gsl_matrix_alloc(params.k_times_d, params.k_times_d_times_s);
  params.st   = gsl_matrix_alloc(params.m_times_d, params.n_times_d);


  
  




  /* vectorize x row-wise */
  x_vec = gsl_vector_view_array(x->data, params.n_times_d);

  /* initiaalize the optimization method */
  switch (method) {
  case 'l': /* LM */
    fdflm.f      = &stls_f_new;
    fdflm.df     = &stls_df_new;
    fdflm.fdf    = &stls_fdf_new;
    fdflm.n      = params.m_times_d;
    fdflm.p      = params.n_times_d;
    fdflm.params = &params;
    
    solverlm = gsl_multifit_fdfsolver_alloc
      (Tlm, params.m_times_d, params.n_times_d);
    gsl_multifit_fdfsolver_set(solverlm, &fdflm, &x_vec.vector);
    g = gsl_vector_alloc(params.n_times_d);
    break;
  case 'q': /* QN 
    fdfqn.f      = &stls_f_;
    fdfqn.df     = &stls_df_;
    fdfqn.fdf    = &stls_fdf_;
    fdfqn.n      = params.n_times_d;
    fdfqn.params = &params;

    solverqn = gsl_multimin_fdfminimizer_alloc(Tqn, params.n_times_d);
    gsl_multimin_fdfminimizer_set(solverqn, &fdfqn, &x_vec.vector, stepqn, opt->epsabs);
    status_dx = GSL_CONTINUE; */
    break;
  case 'n': /* NM */
    fnm.f = &stls_f_;
    fnm.n = params.n_times_d;
    fnm.params = &params;

    solvernm = gsl_multimin_fminimizer_alloc( Tnm, params.n_times_d );
    stepnm = gsl_vector_alloc( params.n_times_d );
    gsl_vector_set_all( stepnm, 0.001 ); /* ??? */
    gsl_multimin_fminimizer_set( solvernm, &fnm, &x_vec.vector, stepnm );
    status_dx   = GSL_CONTINUE;
    status_grad = GSL_CONTINUE;
    break;
  default:
    ; /* error */
  }

  /* optimization loop */
  opt->iter = 0;
  if (opt->disp == 2 || opt->disp == 3)
    PRINTF("STLS optimization:\n");
  do {
    /* print_vec(solverlm->x); */
    opt->iter++;
    switch (method) {
    case 'l': /* Levenberge-Marquardt */
      status = gsl_multifit_fdfsolver_iterate( solverlm );
      /* check for convergence problems */
      if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
	break; /* <- THIS IS WRONG */
      /* check the convergence criteria */
      status_dx = gsl_multifit_test_delta
        (solverlm->dx, solverlm->x, opt->epsabs, opt->epsrel);
      gsl_multifit_gradient(solverlm->J, solverlm->f, g);
      status_grad = gsl_multifit_test_gradient(g, opt->epsgrad);
      /* print information */
      if (opt->disp == 3) {
	gsl_blas_ddot(solverlm->f, solverlm->f, &opt->fmin);
	g_norm = gsl_blas_dnrm2(g);
	PRINTF("%3u: f0 = %16.8f,  ||f0'|| = %16.8f\n", 
	       opt->iter, opt->fmin, g_norm);
      }
      break;
    case 'q':
      status = gsl_multimin_fdfminimizer_iterate( solverqn );
      /* check the convergence criteria */
      status_grad = gsl_multimin_test_gradient(
	   gsl_multimin_fdfminimizer_gradient( solverqn), opt->epsgrad );
      break;
    case 'n':
      status = gsl_multimin_fminimizer_iterate( solvernm );
      /* check the convergence criteria */
      size = gsl_multimin_fminimizer_size( solvernm );
      gsl_multimin_test_size( size, opt->epsabs );
      /* print information */
      if (opt->disp == 3) {
	opt->fmin = gsl_multimin_fminimizer_minimum( solvernm );
	PRINTF("%3u: f0 = %16.8f\n", opt->iter, opt->fmin);
      }
      break;
    }
  } while (status_dx == GSL_CONTINUE && status_grad == GSL_CONTINUE && 
           opt->iter < opt->maxiter);

  switch (method) {
  case  'l':
  /* return the results */
  gsl_vector_memcpy(&x_vec.vector, solverlm->x);
  gsl_multifit_covar(solverlm->J, opt->epsrel, v);
  /* assign the opt output fields */
  gsl_blas_ddot(solverlm->f, solverlm->f, &opt->fmin);
  opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;
  if (opt->iter >= opt->maxiter)
    status = EITER;
  break;
  case 'q':
    break;
  case 'n':
  /* return the results */
  gsl_vector_memcpy(&x_vec.vector, gsl_multimin_fminimizer_x( solvernm ));
  /* gsl_multifit_covar( J??, opt->epsrel, v); */
  /* assign the opt output fields */
  opt->fmin = gsl_multimin_fminimizer_minimum( solvernm );
  opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;
  if (opt->iter >= opt->maxiter)
    status = EITER;
  break;
  }

  /* print exit information */  
  if (opt->disp != 4) { /* unless "off" */
    switch (status) {
    case EITER: 
      PRINTF("STLS optimization terminated by reaching the maximum number " 
	     "of iterations.\nThe result could be far from optimal.\n");
      break;
    case GSL_ETOLF:
      PRINTF("Lack of convergence: progress in function value < EPS.\n");
      break;
    case GSL_ETOLX:
      PRINTF("Lack of convergence: change in parameters < EPS.\n");
      break;
    case GSL_ETOLG:
      PRINTF("Lack of convergence: change in gradient < EPS.\n");
      break;
    default:
      break;
    }
    if ( !status && (opt->disp == 2 || opt->disp == 3) ) /* no error and ( final or iter ) */
      if (status_grad == GSL_CONTINUE)
	PRINTF("Optimization terminated by reaching the convergence " 
	       "tolerance for X.\n");
      else if (status_dx == GSL_CONTINUE)
	PRINTF("Optimization terminated by reaching the convergence " 
	       "tolerance for the gradient.\n");
      else
	PRINTF("Optimization terminated by reaching the convergence " 
	       "tolerance for both X and the gradient.\n");
  }

  switch (method) {
  case 'l': /* LM */
    gsl_multifit_fdfsolver_free(solverlm);
    gsl_vector_free(g);
    break;
  case 'q': /* QN */
    gsl_multimin_fdfminimizer_free(solverqn);
    break;
  case 'n': /* NM */
    gsl_multimin_fminimizer_free(solvernm);
    gsl_vector_free(stepnm);
    break;
  }

  /* free the allocated memory for w */
  for (k = 0; k < w.s; k++) 
    gsl_matrix_free(w.a[k]);
  free(w.a);
  
  /* free preallocated memory for computations */
  gsl_matrix_free(params.tmp);
  gsl_matrix_free(params.gamma);
  free(params.dwork);
  free(params.gamma_vec);

  gsl_matrix_free(params.bx_ext);
  gsl_matrix_free(params.brg_tmp);
  gsl_matrix_free(params.brg_gamma);
  free(params.brg_dwork);
  free(params.brg_gamma_vec);
  gsl_vector_free(params.brg_yr);
  gsl_vector_free(params.brg_f);
  free(params.rb2);


  gsl_vector_free(params.yr);
  free(params.rb);
  gsl_matrix_free(params.x_ext);

  free(params.jres1);
  free(params.jres2);
  gsl_matrix_free(params.dgamma);
  gsl_matrix_free(params.st);


  return GSL_SUCCESS; /* <- correct with status */
}

#define P ((stls_opt_data*) params)

#define M (P->m)
#define N (P->n)
#define D (P->d)

/*
#define M (P->a->size1)
#define N (P->a->size2)
#define D (P->b->size2)*/
#define S (P->w->s)
#define SIZE_W (P->w->a[0]->size1)


/* Compute x_ext into params */
static void compute_xext( const gsl_vector* x, void* params ) {
  /* reshape x as an nxd matrix x_mat */
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );

  /* Form x_ext */
  xmat2xext( x_mat, P->x_ext, params );
}

/* Compute bx_ext into params */
static void compute_bxext( const gsl_vector* x, void* params ) {
  /* reshape x as an nxd matrix x_mat */
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );

  /* Form x_ext */
  xmat2bxext( x_mat, P->bx_ext, params );
}


/* compute f = vec((ax-b)') */
static void compute_f( gsl_vector* f, const gsl_vector* x, void* params ) {

  gsl_matrix_view f_mat = gsl_matrix_view_vector(f, M, D); 
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );
 
  gsl_matrix_memcpy(&f_mat.matrix, P->b);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		 P->a, &x_mat.matrix, -1.0, &f_mat.matrix);
}

static void compute_c_minus_1_2_f( gsl_vector* f, int trans, void* params ) {
  int info;
 /* compute the cost function from the Choleski factor */
  dtbtrs_("U", (trans ? "T" : "N"), "N", 
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &P->one, 
	  P->rb, 
	  &P->k_times_d_times_s, 
	  f->data, 
	  &P->m_times_d, 
	  &info); 
}




static void compute_c_minus_1_f( gsl_vector* f, void* params ) {
  int info;
  
  /* solve for yr by forward substitution */
  dpbtrs_("U",
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &P->one,
	  P->rb,
	  &P->k_times_d_times_s, 
	  f->data, 
	  &P->m_times_d, 
	  &info);
}

static void compute_reshaped_c_minus_1_2_f( gsl_vector* f, int trans, void* params ) {
  int info;
  int i;

  /* compute the cost function from the Choleski factor for each big block of reshaped matrix */
  for (i = 0; i < P->k; i++) {
    dtbtrs_("U", (trans ? "T" : "N"), "N", 
	    &P->d_times_m_div_k, 
	    &P->d_times_s_minus_1, 
	    &P->one, 
	    P->brg_rb, 
	    &P->d_times_s, 
	    (f->data +  i * P->d_times_m_div_k), 
	    &P->d_times_m_div_k, 
	    &info);
	}    
}

static void compute_reshaped_c_minus_1_f( gsl_vector* f, void* params ) {
  int info;
  int i;
  
  /* solve for yr by forward substitution */
  for (i = 0; i < P->k; i++) {
    dpbtrs_("U", 
	    &P->d_times_m_div_k, 
	    &P->d_times_s_minus_1, 
	    &P->one, 
	    P->brg_rb, 
	    &P->d_times_s, 
	    (f->data +  i * P->d_times_m_div_k), 
	    &P->d_times_m_div_k, 
	    &info);
	}    
}




/*************************************
 * This function convert between  representations
 *     D * K * m_div_k  array (original)
 *     D * m_div_k * K  array (reshaped)
 *
 *  forward = 1  :  original -> reshaped
 *  forward = 0  :  reshaped -> original
 * 
 * 
 * Reshaped vector can be multiplied by (smaller) block of reshaped gamma matrix
 *************************************/
void reshape_f( gsl_vector *reshaped, gsl_vector * original, void* params, int forward ) {
  gsl_vector_view w_orig, w_resh;
   int j, i; 
    
  for (j = 0; j < P->m_div_k;  j++)  {
    for (i = 0; i < P->k;  i++)  {
      w_orig = gsl_vector_subvector(original, (j*P->k + i) *D, D);
      w_resh = gsl_vector_subvector(reshaped, (j + i * P->m_div_k) *D, D);
       
      if (forward)  {
        gsl_vector_memcpy(&w_resh.vector, &w_orig.vector);
      } else {
        gsl_vector_memcpy(&w_orig.vector, &w_resh.vector);
      }
    }
  }  
}








/* 
*  STLS_F: STLS cost function evaluation 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*  f      - cost function
*/

int stls_f (const gsl_vector* x, void* params, gsl_vector* f)
{
  compute_f(f, x, params);
  compute_xext(x, params);
  cholgam(params);
  compute_c_minus_1_2_f(f, 1, params);
 
  return GSL_SUCCESS;
}


int stls_f_new (const gsl_vector* x, void* params, gsl_vector* f)
{
  compute_f(f, x, params);
  compute_bxext(x, params);
  cholbrg(params);

  reshape_f(P->brg_f, f, params, 1);
  compute_reshaped_c_minus_1_2_f(P->brg_f, 1, params);
  reshape_f(P->brg_f, f, params, 0);

  return GSL_SUCCESS;
}



/* 
*  STLS_F_: STLS cost function evaluation for QN 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*/

double stls_f_ (const gsl_vector* x, void* params)
{
  double ftf;

  /* Use yr as a temporary variable */
  stls_f(x, params, P->yr);
  gsl_blas_ddot(P->yr, P->yr, &ftf);

  return ftf;
}


/* 
*  STLS_FD: STLS first derivative evaluation 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*  deriv  - first derivative
*/

int stls_df (const gsl_vector* x, 
	     void* params, gsl_matrix* deriv)
{
  compute_xext(x, params);
  cholgam(params);
  compute_f(P->yr, x, params);

  compute_c_minus_1_f(P->yr, params);
  jacobian(params, deriv);

  return GSL_SUCCESS;
}



int stls_df_new (const gsl_vector* x, void* params, gsl_matrix* deriv)
{
  compute_bxext(x, params);
  cholbrg(params);
  compute_f(P->yr, x, params);
  

  reshape_f(P->brg_yr, P->yr, params, 1);
  compute_reshaped_c_minus_1_f(P->brg_yr, params);
 /* reshape_f(P->brg_yr, P->yr, params, 0);


  cholgam(params); */

  jacobian_new(params, deriv);

  return GSL_SUCCESS;
}




/* 
*  STLS_FDF: STLS cost function and first derivative evaluation
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*  f      - cost function
*  deriv  - first derivative
*/

int stls_fdf (const gsl_vector* x, void* params, gsl_vector* f, 
	      gsl_matrix* deriv)
{
  compute_f(f, x, params);
  compute_xext(x, params);
  cholgam(params);
  compute_c_minus_1_2_f(f, 1, params);
  gsl_vector_memcpy(P->yr, f);
  compute_c_minus_1_2_f(P->yr, 0, params);
  jacobian(params, deriv);

  return GSL_SUCCESS;
}

int stls_fdf_new (const gsl_vector* x, void* params, gsl_vector* f, 
	      gsl_matrix* deriv)
{
  compute_f(P->yr, x, params);
  compute_bxext(x, params);
  cholbrg(params);

  reshape_f(P->brg_yr, P->yr, params, 1);
  compute_reshaped_c_minus_1_2_f(P->brg_yr, 1, params);
  reshape_f(P->brg_yr, f, params, 0);

  compute_reshaped_c_minus_1_2_f(P->brg_yr, 0, params);
  
  jacobian_new(params, deriv);

  return GSL_SUCCESS;
}




/* ---- Auxiliary functions ---- */

/* 
*  XMAT2XEXT: x_mat |-> x_ext
*
*  x_mat  - the parameter x viewed as a matrix, 
*  x_ext  = kron( I_k, [ x; -I_d ] ),
*  params - used for the dimensions n = rowdim(X), 
*           d = coldim(X), k, and n_plus_d = n + d
*/



void xmat2bxext( gsl_matrix_const_view x_mat, gsl_matrix *bx_ext,  stls_opt_data* params )
{
  gsl_vector_view diag;
  gsl_matrix_view submat;

  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  submat = gsl_matrix_submatrix(bx_ext, 0, 0, N, D); /* select x in (1,1) */
  gsl_matrix_memcpy(&submat.matrix, &x_mat.matrix); /* assign x */
  submat = gsl_matrix_submatrix(bx_ext, N, 0, D, D); /* select -I in (1,1)*/
  diag   = gsl_matrix_diagonal(&submat.matrix);     /* assign -I */
  gsl_vector_set_all(&diag.vector, -1);
}



void xmat2xext( gsl_matrix_const_view x_mat, gsl_matrix *x_ext,  stls_opt_data* params )
{
  int i;
  gsl_matrix_view submat, source;

  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  xmat2bxext(x_mat, x_ext, params);

  /* copy block (1,1) in (2,2),...,(k,k) */
  source = gsl_matrix_submatrix(x_ext, 0, 0, P->n_plus_d, D);
  for (i = 1; i < P->k; i++) {
    submat = gsl_matrix_submatrix(x_ext, i*P->n_plus_d, i*D, P->n_plus_d, D);
    gsl_matrix_memcpy(&submat.matrix, &source.matrix);
  }
}






/* 
*  CHOLGAM: Choleski factorization of the block of reshaped Gamma  matrix
*
*  params - used for the dimensions n = rowdim(X), 
*           d = coldim(X), k, and n_plus_d = n + d
*
*  params.x_ext  = kron( I_k, [ x; -I_d ] ),
*
*  Output:
*    params.brg_rb     - Cholesky factor in a packed form
*/

void cholbrg( stls_opt_data* params )
{
  int k, info;
  gsl_matrix_view submat;
  gsl_matrix_view  b_w_k; 
  gsl_matrix *bx_ext = P->bx_ext;



  const int zero = 0;
  /* MB02GD has very bad description of parameters */


  /* Take block of x_ext - b_x_ext = x_ext(1:(n+d),1:d) * /
  gsl_matrix_view b_x_ext;
  b_x_ext =  gsl_matrix_submatrix(P->x_ext, 0, 0, P->n_plus_d, D);
  bx_ext = &b_x_ext.matrix; */
  

  /* compute brgamma_k = b_x_ext' * w_k(1:(n+d),1:(n+d)) * b_x_ext */
  for (k = 0; k < S; k++) {
    submat = gsl_matrix_submatrix(P->brg_gamma, 0, k*D, D, D);
    b_w_k = gsl_matrix_submatrix(P->w->a[k], 0, 0, P->n_plus_d, P->n_plus_d);
				  
				  
    /* compute tmp = x_ext' * w_k */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, bx_ext, &b_w_k.matrix, 0.0, P->brg_tmp);
    /* compute brgamma_k = tmp * x_ext */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,  P->brg_tmp, bx_ext, 0.0, &submat.matrix);
  }

  gsl_matrix_vectorize(P->brg_gamma_vec, P->brg_gamma);
  
/*  PRINTF("cholgam: PDW = %d, LDWORK = %d, 3K = %d\n, s-1 = %d", 1 + (P->s_minus_1 + 1)* P->k_times_d * P->k_times_d, ldwork, 3 * P->k_times_d, P->s_minus_1 );*/
  /* Cholesky factorization of Gamma */
  mb02gd_("R", "N",
	  &P->d, 	/* block size */
	  &P->m_div_k,		/* block_dim(brGAMMA) */
	  &P->s_minus_1, 	/* non-zero block superdiagonals */
	  &zero, 		/* previous computations */
	  &P->m_div_k, 	        /* to-be-computed */
	  P->brg_gamma_vec, /* non-zero part of first block row of brGAMMA */
	  &P->d,      	/* row_dim(brg_gamma) */
	  P->brg_rb, 			/* packed Cholesky factor */
	  &P->d_times_s, /* row_dim(rb) */
	  P->brg_dwork, &P->brg_ldwork, &info); /**/


  /* check for errors of mb02gd */
  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED */
  }
}


/* 
*  CHOLGAM: Choleski factorization of the Gamma matrix
*
*  params - used for the dimensions n = rowdim(X), 
*           d = coldim(X), k, and n_plus_d = n + d
*
*  params.x_ext  = kron( I_k, [ x; -I_d ] ),
*
*  Output:
*    params.rb     - Cholesky factor in a packed form
*/

void cholgam( stls_opt_data* params )
{
  int k, info;
  gsl_matrix_view submat, source;



  const int zero = 0;
  /* MB02GD has very bad description of parameters */

  /* compute gamma_k = x_ext' * w_k * x_ext */
  for (k = 0; k < S; k++) {
    submat = gsl_matrix_submatrix(P->gamma, 0, 
				  k*P->k_times_d, P->k_times_d, P->k_times_d);
    /* compute tmp = x_ext' * w_k */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, 
		   P->x_ext, P->w->a[k], 0.0, P->tmp);
    /* compute tmp * x_ext */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		   P->tmp, P->x_ext, 0.0, &submat.matrix);
  }

  gsl_matrix_vectorize(P->gamma_vec, P->gamma);
  
  /*PRINTF("cholgam: PDW = %d, LDWORK = %d, 3K = %d\n, s-1 = %d", 1 + (P->s_minus_1 + 1)* P->k_times_d * P->k_times_d, ldwork, 3 * P->k_times_d, P->s_minus_1 );*/

  /* Cholesky factorization of Gamma */
  mb02gd_("R", "N",  
	  &P->k_times_d, 	/* block size */
	  &P->m_div_k,		/* block_dim(GAMMA) */
	  &P->s_minus_1, 	/* non-zero block superdiagonals */
	  &zero, 		/* previous computations */
	  &P->m_div_k, 	        /* to-be-computed */
	  P->gamma_vec, /* non-zero part of first block row of GAMMA */
	  &P->k_times_d,      	/* row_dim(gamma) */
	  P->rb, 			/* packed Cholesky factor */
	  &P->k_times_d_times_s, /* row_dim(rb) */
	  P->dwork, &P->ldwork, &info);


  /* check for errors of mb02gd */
  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED */
  }
}


/* 
*  JACOBIAN: Computes the Jacobian 
*
*  params - used for the dimensions n = rowdim(X), 
*    .x_ext  = kron( I_k, [ x; -I_d ] ),
*           d = coldim(X), k, and n_plus_d = n + d
*    .rb     - Choleski factor in a packed form
*    .yr     - solution of the system Gamma*yr = r
* output
*  deriv  - Jacobian
*/

void jacobian( stls_opt_data* params, gsl_matrix* deriv )
{
  int i, j, l, k, info;
  const int zero = 0, one = 1;
  gsl_matrix_view submat; 
  gsl_vector_view subvec, res_vec;

	gsl_vector_view tmp1_row, tmp1_col;
	gsl_vector_view w_k_row, w_k_col;


  /* first term of the Jacobian Gamma^{-1/2} kron(a,I_d) */

  /* deriv = kron(a,I_d) */
  gsl_matrix_set_zero(deriv);
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
      submat = gsl_matrix_submatrix(deriv, i*D, j*D, D, D);
      /* assign a(i,j) on the diagonal of submat */
      for (l = 0; l < D; l++)
				gsl_matrix_set(&submat.matrix, l, l, gsl_matrix_get(P->a, i, j));
    }
  
  /* res1 = vec(deriv), FORTRAN column-major convention */
  gsl_matrix_vectorize(P->jres1, deriv);

  /* res =  Gamma^{-1/2} * res  */
  dtbtrs_("U", "T", "N", 
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &deriv->size2, 
	  P->rb, 
	  &P->k_times_d_times_s, 
	  P->jres1, 
	  &P->m_times_d, 
	  &info);
  
  /* convert back deriv = res1 in the C row-major convention */
  gsl_matrix_vec_inv(deriv, P->jres1);


  /* second term (naive implementation) */
  res_vec = gsl_vector_view_array(P->jres2, P->m_times_d);

  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
			/* form dgamma = d/d x_ij gamma */
			gsl_matrix_set_zero(P->dgamma);
			for (k = 0; k < S; k++) {
	  	  submat = gsl_matrix_submatrix(P->dgamma, 0, k*P->k_times_d, P->k_times_d, P->k_times_d);

  	  	/* compute tmp1 = dx_ext' * w_k * x_ext * /
		    /* Iterate over rows of dx_ext' * w_k */
			  for (l = 0; l < P->k; l++) {
			    tmp1_row = gsl_matrix_row (&submat.matrix, l*D+j);
			    w_k_row = gsl_matrix_row (P->w->a[k], l*P->n_plus_d+i);
			    gsl_blas_dgemv(CblasTrans, 1.0, P->x_ext, &w_k_row.vector, 0.0, &tmp1_row.vector); 
			  }
			  
    	  /* compute submat = submat  + x_ext' * tmp * dx_ext * /
		    /* Iterate over rows of dx_ext' * w_k' */
			  for (l = 0; l < P->k; l++) {
			    tmp1_col = gsl_matrix_column (&submat.matrix, l*D+j);
			    w_k_col = gsl_matrix_column (P->w->a[k], l*P->n_plus_d+i);
			    gsl_blas_dgemv(CblasTrans, 1.0, P->x_ext, &w_k_col.vector, 1.0, &tmp1_col.vector); 
			  }
     	}
    	/* compute st_ij = DGamma * yr */
    	subvec = gsl_matrix_column(P->st, i*D+j);
    	tmv_prod(P->dgamma, S, P->yr, P->m_div_k, &subvec.vector);
    	gsl_vector_memcpy(&res_vec.vector, &subvec.vector);
    	/* solve st_ij = Gamma^{-1/2}st_ij */
    	dtbtrs_("U", "T", "N", 
					&P->m_times_d, 
					&P->k_times_d_times_s_minus_1, 
					&one, 
					P->rb, 
					&P->k_times_d_times_s, 
					P->jres2,
					&P->m_times_d, 
				&info);
				gsl_vector_memcpy(&subvec.vector, &res_vec.vector);
    }
  }




  /* deriv = deriv - 0.5 * st */
  gsl_matrix_scale(P->st, 0.5);
  gsl_matrix_sub(deriv, P->st);

}




void jacobian_new( stls_opt_data* params, gsl_matrix* deriv )
{
  int i, j, l, k, info;
  const int zero = 0, one = 1;
  gsl_matrix_view submat, mat; 
  gsl_vector_view subvec, res_vec;

    
  
  
  

	gsl_vector_view tmp1_row, tmp1_col;
	gsl_vector_view w_k_row, w_k_col;


  /* first term of the Jacobian Gamma^{-1/2} kron(a,I_d) */

  /* deriv = kron(a,I_d) * /
  gsl_matrix_set_zero(deriv);
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
      submat = gsl_matrix_submatrix(deriv, i*D, j*D, D, D);
      /* assign a(i,j) on the diagonal of submat * /
      for (l = 0; l < D; l++)
				gsl_matrix_set(&submat.matrix, l, l, gsl_matrix_get(P->a, i, j));
    }
    
  
  /* res1 = vec(deriv), FORTRAN column-major convention * /
  gsl_matrix_vectorize(P->jres1, deriv);

  /* res =  Gamma^{-1/2} * res  * /
  dtbtrs_("U", "T", "N", 
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &deriv->size2, 
	  P->rb, 
	  &P->k_times_d_times_s, 
	  P->jres1, 
	  &P->m_times_d, 
	  &info);
  
  /* convert back deriv = res1 in the C row-major convention * /
  gsl_matrix_vec_inv(deriv, P->jres1); */

  double *test = calloc(P->d_times_m_div_k * D, sizeof(double));
  gsl_matrix *mattemp = gsl_matrix_calloc(P->d_times_m_div_k, D);

  gsl_matrix_view submat1, submat2;
  int ibr, ik, i1, j1, kbr;
  double a_kj;
  
  gsl_matrix_set_zero(P->st);
  
  for (j = 0; j < N; j++) {
    for (ik = 0; ik < P->k; ik++) {

      /* Fill right-hand matrix */
      for (ibr = 0; ibr < P->m_div_k; ibr++) {
        for (k = 0; k < D; k++) {
          for (l = 0; l < D ;l++) {
            gsl_matrix_set(mattemp, k + ibr*D, l, 0);
          }
          
          
          gsl_matrix_set(mattemp, k + ibr*D, k, gsl_matrix_get(P->a, ik + ibr * P->k,j)) ;
        }
      }
      

      gsl_matrix_vectorize(test, mattemp);
      
      
      dtbtrs_("U", "T", "N", 
	      &P->d_times_m_div_k, 
	      &P->d_times_s_minus_1, 
	      &P->d, 
  	    P->brg_rb, 
	      &P->d_times_s, 
	      test, 
	      &P->d_times_m_div_k, 
  	    &info);

      gsl_matrix_vec_inv(mattemp, test);


      for (ibr = 0; ibr < P->m_div_k; ibr++) {
        submat1 = gsl_matrix_submatrix(P->st, (ik + ibr*P->k) * D, j*D, D, D); 
        submat2 = gsl_matrix_submatrix(mattemp, ibr*D, 0, D, D); 
        
        gsl_matrix_memcpy(&submat1.matrix, &submat2.matrix);
			}
    }
  }
    
  free(test);
  gsl_matrix_free(mattemp);
  
  

  gsl_matrix_memcpy(deriv, P->st);
/*  PRINTF("jacobian_new: step1 max, min = %f, %f, \n", gsl_matrix_max(P->st), gsl_matrix_min(P->st));*/





  /* second term (naive implementation) */
  res_vec = gsl_vector_view_array(P->jres2, P->m_times_d);

  gsl_matrix *bx_ext = P->bx_ext;


  /* Take block of x_ext - b_x_ext = x_ext(1:(n+d),1:d) * /
  gsl_matrix_view b_x_ext;
  b_x_ext =  gsl_matrix_submatrix(P->x_ext, 0, 0, P->n_plus_d, D);
  bx_ext = &b_x_ext.matrix;  */


  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {

      /* Alternative */

		  /* form dgamma = d/d x_ij gamma */
			gsl_matrix_set_zero(P->dgamma);


			for (k = 0; k < S; k++) {
			  gsl_matrix_view  b_w_k = gsl_matrix_submatrix(P->w->a[k], 0, 0, P->n_plus_d, P->n_plus_d);
			
	  	  submat = gsl_matrix_submatrix(P->dgamma, 0, k*P->d, P->d, P->d);

  	  	/* compute tmp1 = dx_ext' * w_k * x_ext * /
		    /* Iterate over rows of dx_ext' * w_k */
		    tmp1_row = gsl_matrix_row (&submat.matrix, j);
		    w_k_row = gsl_matrix_row (&b_w_k.matrix, i);
		    gsl_blas_dgemv(CblasTrans, 1.0, bx_ext, &w_k_row.vector, 0.0, &tmp1_row.vector); 
			  
    	  /* compute submat = submat  + x_ext' * tmp * dx_ext * /
		    /* Iterate over rows of dx_ext' * w_k' */
  	    tmp1_col = gsl_matrix_column (&submat.matrix, j);
		    w_k_col = gsl_matrix_column (&b_w_k.matrix, i);
		    gsl_blas_dgemv(CblasTrans, 1.0, bx_ext, &w_k_col.vector, 1.0, &tmp1_col.vector); 
     	}
     	
     	submat = gsl_matrix_submatrix(P->dgamma, 0, 0, P->d, S*P->d);
     	
    	/* compute st_ij = DGamma * yr */
  
      subvec = gsl_matrix_column(P->st, i*D+j);
    		
      for (l = 0; l < P->k; l++) { 
        gsl_vector_view subvec_res = gsl_vector_subvector(&res_vec.vector, l * P->d_times_m_div_k, P->d_times_m_div_k);
        gsl_vector_view subvec_yr = gsl_vector_subvector(P->brg_yr, l * P->d_times_m_div_k, P->d_times_m_div_k);
				
				
        tmv_prod(&submat.matrix, S, &subvec_yr.vector, P->m_div_k, &subvec_res.vector);  

      	/* solve st_ij = Gamma^{-1/2}st_ij */
      	dtbtrs_("U", "T", "N", 
			  		&P->d_times_m_div_k, 
				  	&P->d_times_s_minus_1, 
  					&one, 
	  				P->brg_rb, 
		  			&P->d_times_s, 
			  		subvec_res.vector.data,
				  	&P->d_times_m_div_k, 
    				&info);
    				
    				
    	/*	double sum_max = 0, sum_min = 0;		*/
    		for (k = 0; k < P->m_div_k; k++) {
    		  gsl_vector_view subvec1 = gsl_vector_subvector(&subvec_res.vector, D * k, D);
    		  gsl_vector_view subvec2 = gsl_vector_subvector(&subvec.vector, k * P->k_times_d + l * D, D);
    		 
    		  gsl_vector_memcpy(&subvec2.vector, &subvec1.vector);
    		 
    		 /* gsl_vector_sub(&subvec1.vector, &subvec2.vector);
    		  gsl_vector_mul(&subvec1.vector, &subvec1.vector);
    		  
    		  sum_max += gsl_vector_max(&subvec1.vector);*/
    		}
    		
/*    		PRINTF("jacobian_new: max = %14.10f\n", sum_max);*/
    	}

			
			
    }
  }




  /* deriv = deriv - 0.5 * st */
  gsl_matrix_scale(P->st, 0.5);
  gsl_matrix_sub(deriv, P->st);

}



/*
* tmv_prod: block-Toeplitz matrix T times vector v
*
* t - nonzero part of the first block row of T
* s - number of blocks in t, t = [t_0 ... t_s-1]
* s_1 = s - 1;  m = (int) v->size1 / t->size1
* p - result
*/ 

void tmv_prod(gsl_matrix* t, int s, gsl_vector* v, int m, 
	      gsl_vector* p)
{
  int i, imax, temp, s_1 = s - 1;
  gsl_vector_view subv, subp; 	/* subvector of v and p */
  gsl_matrix* tt;		/* = [t_s-1' ... t_1' t_0 t_1 ... t_s-1] */
  gsl_matrix_view submat, source;

#define TM (t->size1) 		/* = block size */
#define TN (t->size2)		/* = s(block size) */

  /* form tt */
  tt = gsl_matrix_alloc(TM, 2*TN - TM);
  for (i = 0; i < s_1; i++) {
    submat = gsl_matrix_submatrix(tt, 0, i*TM, TM, TM);
    source = gsl_matrix_submatrix(t, 0, (s_1-i)*TM, TM, TM);
    gsl_matrix_transpose_memcpy
      (&submat.matrix, &source.matrix);
  }
  submat = gsl_matrix_submatrix(tt, 0, s_1*TM, TM, TN);
  gsl_matrix_memcpy(&submat.matrix, t);

  /* construct p = T*v */
  gsl_vector_set_zero(p);

  /* beginning and end parts of the product p */
  for (i = 0; i < s_1; i++) {
    temp = (s+i)*TM;
    /* beginning part */
    subp = gsl_vector_subvector(p, i*TM, TM);
    subv = gsl_vector_subvector(v, 0, temp);
    submat = gsl_matrix_submatrix
      (tt, 0, (s_1-i)*TM, TM, temp);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &submat.matrix, 
		   &subv.vector, 0.0, &subp.vector);
    /* last part */
    subp = gsl_vector_subvector(p, p->size - (i+1)*TM, TM);
    subv = gsl_vector_subvector(v, v->size - temp, temp);
    submat = gsl_matrix_submatrix(tt, 0, 0, TM, temp);    
    gsl_blas_dgemv(CblasNoTrans, 1.0, &submat.matrix, 
		   &subv.vector, 0.0, &subp.vector);
  }

  /* middle part */
  for (i = s_1, imax = m - s_1 ; i < imax; i++) {
    subp = gsl_vector_subvector(p, i*TM, TM);
    subv = gsl_vector_subvector(v, (i-s_1)*TM, tt->size2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, tt, &subv.vector, 
		   0.0, &subp.vector);
  }
  
  
  gsl_matrix_free(tt);
}


/* s2w: finds the covariance matrices w from the data structure */

#define NL  (s->a[l].ncol) 	/* # col. in C^l */
#define NBL (s->a[l].nb) 	/* # col. in a repeated block of C^l */

int s2w(const data_struct* s, w_data* w)
{
  int n_d, k_n_d, k, l, i, offset, imax, sum_nl;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;
  char err_msg[70];
  
  /* find n_d = n+d = sum_{l=1}^q n_l and w->s */
  for (n_d = l = 0, w->s = 1; l < s->q; l++) {
    n_d = n_d + NL;
    if ((s->a[l].type == 'T' || s->a[l].type == 'H') && NL/NBL > w->s)
      w->s = NL/NBL;
  }
  k_n_d = s->k * n_d;

  w->a = (gsl_matrix**) malloc(w->s * sizeof(gsl_matrix *));
  zk   = gsl_matrix_alloc(n_d, n_d);
  /* construct w */
  for (k = 0; k < w->s; k++) { 
    gsl_matrix_set_zero(zk);
    for (l = sum_nl = 0; l < s->q; sum_nl += s->a[l++].ncol) { 
      zkl = gsl_matrix_submatrix(zk, sum_nl, sum_nl, NL, NL); 
      switch (s->a[l].type) {
      case 'T': case 'H':
	offset = NBL*k;
	imax   = NL - offset;
	for (i = 0; i < imax; i++)
	  if (s->a[l].type == 'H')
	    gsl_matrix_set(&zkl.matrix, i+offset, i, 1);
	  else
	    gsl_matrix_set(&zkl.matrix, i, i+offset, 1);
	break;
      case 'U': 
	if (k == 0)
	  gsl_matrix_set_identity(&zkl.matrix);
	/* else zik is a zero matrix */
	break;
      case 'E': 
	/* zik is a zero matrix */
	break;
      default:
	sprintf(err_msg, "Unknown structure type %c.",
		s->a[l].type);
	GSL_ERROR(err_msg, GSL_EINVAL);
	break;
      }
    }
    w->a[k] = gsl_matrix_calloc(k_n_d, k_n_d);
    /* w->a[k] = kron(Ik, zk) */
    for (i = 0; i < s->k; i++) {
      /* select the i-th diagonal block in a matrix view */
      wi = gsl_matrix_submatrix( w->a[k], i*n_d, 
				 i*n_d, n_d, n_d );
      gsl_matrix_memcpy(&wi.matrix, zk);
    }
  }
  gsl_matrix_free(zk);

  return GSL_SUCCESS;
}

/* #include "corr.c" TO BE ADDED */

/* ********************************************************************* */
/* Interface to the SLICOT function MB02MD. Solves a TLS problem AX = B. */
/* ********************************************************************* */

int tls(gsl_matrix* a, gsl_matrix* b, gsl_matrix* x)
{
  double *c, *xa, *s, *dwork, tol=0;
  int    m, n, l, iwarn, info, ldwork, *iwork, ldc, r=0; 

  /* Define constants */
  m = a->size1;
  n = a->size2;
  l = b->size2;			/* = d */
  ldc = GSL_MAX(m, n+l);
  if ( m < n+l ) {
    ldwork = GSL_MAX(m*(n+l)+GSL_MAX(3*m+n+l, 5*m), 3*l);
  } else {
    ldwork = GSL_MAX(3*(n+l)+m, 5*(n+l));
  }
  
  /* c = [a b] in FORTRAN format */
  c = malloc( m * (n+l) * sizeof(double));
  gsl_matrix_vectorize(c, a);
  gsl_matrix_vectorize(c + m*n, b);

  xa = malloc( n * l * sizeof(double));
  s  = malloc( (n+l) * sizeof(double));
  dwork = malloc( ldwork * sizeof(double));
  iwork = malloc( l * sizeof(int));

  /* Call MB02MD */
  mb02md_("B", &m, &n, &l, &r, c, &ldc, s, xa, &n, &tol, iwork, dwork, &ldwork, &iwarn, &info);

  /* copy the result in x */
  if (!info)
    gsl_matrix_vec_inv(x, xa);

  /* free memory */
  free(c);
  free(xa);
  free(s);
  free(dwork);
  free(iwork);

  return info;
}

/* gsl_matrix_vectorize: vectorize column-wise a gsl_matrix */
void gsl_matrix_vectorize(double* v, gsl_matrix* m)
{
  int i, j;

  for (i = 0; i < m->size1; i++)
    for (j = 0; j < m->size2; j++)
      v[i+j*m->size1] = gsl_matrix_get(m,i,j);
}


/* gsl_matrix_vec_inv: gsl_matrix from an array */
void gsl_matrix_vec_inv(gsl_matrix* m, double* v)
{
  int i, j;

  for (i = 0; i < m->size1; i++)
    for (j = 0; j < m->size2; j++)
      gsl_matrix_set(m,i,j,v[i+j*m->size1]);
}


/* print matrix */
void print_mat(const gsl_matrix* m)
{
  int i, j;

  PRINTF("\n");
  for (i = 0; i < m->size1; i++) {
    for (j = 0; j < m->size2; j++)
      PRINTF("%16.14f ", gsl_matrix_get(m, i, j));
    PRINTF("\n");
  }
  PRINTF("\n");
}


/* print matrix */
void print_mat_tr(const gsl_matrix* m)
{
  int i, j;

  PRINTF("\n");
  for (j = 0; j < m->size2; j++) {
    for (i = 0; i < m->size1; i++) {
      PRINTF("%16.14f ", gsl_matrix_get(m, i, j));
    }
    PRINTF("\n");
  }
  PRINTF("\n");
}


/* print_arr: print array */
void print_arr(double* a, int n)
{
  int i;

  PRINTF("\n");
  for (i = 0; i < n; i++)
    PRINTF("%f\n",*(a+i));
  PRINTF("\n");
}
