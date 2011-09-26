/* stls.c: implementations of the functions from stls.h */

#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include "stls.h"

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
  params.size_of_gamma = params.k_times_d * 
    params.k_times_d_times_s * sizeof(double);
  params.size_of_rb = params.m_times_d * 
    params.k_times_d_times_s * sizeof(double);
  params.size_of_dwork = LDWORK * sizeof(double);

  /* vectorize x row-wise */
  x_vec = gsl_vector_view_array(x->data, params.n_times_d);

  /* initiaalize the optimization method */
  switch (method) {
  case 'l': /* LM */
    fdflm.f      = &stls_f;
    fdflm.df     = &stls_df;
    fdflm.fdf    = &stls_fdf;
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

  return GSL_SUCCESS; /* <- correct with status */
}

#define P ((stls_opt_data*) params)
#define M (P->a->size1)
#define N (P->a->size2)
#define D (P->b->size2)
#define S (P->w->s)
#define SIZE_W (P->w->a[0]->size1)

/* 
*  STLS_F: STLS cost function evaluation 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*  f      - cost function
*/

int stls_f (const gsl_vector* x, void* params, gsl_vector* f)
{
  int info;
  const int one = 1;
  gsl_matrix *x_ext;
  gsl_matrix_view submat; 
  double *rb;

  /* reshape x as an nxd matrix x_mat */
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );

  /* reserve memory for x_ext and form it */
  x_ext = gsl_matrix_calloc(SIZE_W, P->k_times_d); /* SIZE_W = ndk */
  xmat2xext( x_mat, x_ext, params );

  /* Cholesky factorization */
  rb = (double*) malloc(P->size_of_rb);
  cholgam(x_ext, params, rb);

  /* compute f = vec((ax-b)') */
  submat = gsl_matrix_view_vector(f, M, D);
  gsl_matrix_memcpy(&submat.matrix, P->b);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		 P->a, &x_mat.matrix, -1.0, &submat.matrix);
  
  /* compute the cost function from the Choleski factor */
  dtbtrs_("U", "T", "N", 
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &one, 
	  rb, 
	  &P->k_times_d_times_s, 
	  f->data, 
	  &P->m_times_d, 
	  &info);

  free(rb);
  gsl_matrix_free(x_ext);

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
  gsl_vector* f;
  double ftf;

  f = gsl_vector_alloc( ((stls_opt_data*) params)->m_times_d );
  stls_f(x, params, f);
  gsl_blas_ddot(f, f, &ftf);

  gsl_vector_free(f);
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
  int info;
  const int one = 1;
  double* rb;
  gsl_matrix *x_ext;
  gsl_matrix_view submat; 
  gsl_vector *yr;

  /* reshape x as an nxd matrix x_mat */
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );

  /* reserve memory for x_ext and form it */
  x_ext = gsl_matrix_calloc(SIZE_W, P->k_times_d); /* SIZE_W = ndk */
  xmat2xext( x_mat, x_ext, params );

  /* Cholesky factorization */
  rb = (double*) malloc(P->size_of_rb);
  cholgam(x_ext, params, rb);
  
  /* compute vec((ax-b)') */
  yr = gsl_vector_alloc(P->m_times_d);
  submat = gsl_matrix_view_vector(yr, M, D);
  gsl_matrix_memcpy(&submat.matrix, P->b);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		 P->a, &x_mat.matrix, -1.0, &submat.matrix);

  /* solve for yr by forward substitution */
  dpbtrs_("U",
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &one,
	  rb,
	  &P->k_times_d_times_s, 
	  yr->data, 
	  &P->m_times_d, 
	  &info);

  jacobian(x_ext, params, rb, yr, deriv);

  free(rb);

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
  int info;
  const int one = 1;
  double* rb;
  gsl_matrix *x_ext;
  gsl_matrix_view submat; 
  gsl_vector *yr;

  /* reshape x as an nxd matrix x_mat */
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D);     

  /* reserve memory for x_ext and form it */
  x_ext = gsl_matrix_calloc(SIZE_W, P->k_times_d); /* SIZE_W = ndk */
  xmat2xext( x_mat, x_ext, params );

  /* Cholesky factorization */
  rb = (double*) malloc(P->size_of_rb);
  cholgam(x_ext, params, rb);

  /* compute f = vec((ax-b)') */
  submat = gsl_matrix_view_vector(f, M, D);
  gsl_matrix_memcpy(&submat.matrix, P->b);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		 P->a, &x_mat.matrix, -1.0, &submat.matrix);

  /* compute the cost function from the Choleski factor */
  dtbtrs_("U", "T", "N", 
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &one, 
	  rb, 
	  &P->k_times_d_times_s, 
	  f->data, 
	  &P->m_times_d, 
	  &info);

  /* compute vec((ax-b)') */
  yr = gsl_vector_alloc(P->m_times_d);
  gsl_vector_memcpy(yr, f);

  /* solve for yr by forward substitution */
  dtbtrs_("U", "N", "N", 
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &one, 
	  rb,
	  &P->k_times_d_times_s, 
	  yr->data, 
	  &P->m_times_d, 
	  &info);

  jacobian( x_ext, params, rb, yr, deriv );

  free(rb);

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

void xmat2xext( gsl_matrix_const_view x_mat,
	        gsl_matrix* x_ext, 
	        stls_opt_data* params )
{
  int i;
  gsl_vector_view diag;
  gsl_matrix_view submat, source;

  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  submat = gsl_matrix_submatrix(x_ext, 0, 0, N, D); /* select x in (1,1) */
  gsl_matrix_memcpy(&submat.matrix, &x_mat.matrix); /* assign x */
  submat = gsl_matrix_submatrix(x_ext, N, 0, D, D); /* select -I in (1,1)*/
  diag   = gsl_matrix_diagonal(&submat.matrix);     /* assign -I */
  gsl_vector_set_all(&diag.vector, -1);

  /* copy block (1,1) in (2,2),...,(k,k) */
  source = gsl_matrix_submatrix(x_ext, 0, 0, P->n_plus_d, D);
  for (i = 1; i < P->k; i++) {
    submat = gsl_matrix_submatrix(x_ext, i*P->n_plus_d, i*D, P->n_plus_d, D);
    gsl_matrix_memcpy(&submat.matrix, &source.matrix);
  }
}


#define mymax(a,b) ((a) > (b) ? (a) : (b)) 

/* 
*  CHOLGAM: Choleski factorization of the Gamma matrix
*
*  x_ext  = kron( I_k, [ x; -I_d ] ),
*  params - used for the dimensions n = rowdim(X), 
*           d = coldim(X), k, and n_plus_d = n + d
*  rb     - Choleski factor in a packed form
*/

void cholgam( gsl_matrix* x_ext, 
	      stls_opt_data* params, double* rb )
{
  int k, info;
  gsl_matrix *tmp, *gamma;
  gsl_matrix_view submat, source;
  double *gamma_vec, *dwork;
  const int zero = 0;
  /* MB02GD has very bad description of parameters */
  const int ldwork = 1 + (P->s_minus_1 + 1)* P->k_times_d * P->k_times_d +  /* pDW */ 
                       3 * P->k_times_d + /* 3 * K */
                       mymax(P->s_minus_1 + 1,  P->m_div_k - 1 - P->s_minus_1) * P->k_times_d * P->k_times_d; /* Space needed for MB02CV */

  /* compute gamma_k = x_ext' * w_k * x_ext */
  gamma = gsl_matrix_alloc(P->k_times_d, P->k_times_d_times_s);
  tmp   = gsl_matrix_alloc(P->k_times_d, SIZE_W);
  for (k = 0; k < S; k++) {
    submat = gsl_matrix_submatrix(gamma, 0, 
				  k*P->k_times_d, P->k_times_d, P->k_times_d);
    /* compute tmp = x_ext' * w_k */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, 
		   x_ext, P->w->a[k], 0.0, tmp);
    /* compute tmp * x_ext */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		   tmp, x_ext, 0.0, &submat.matrix);
  }
  gsl_matrix_free(tmp);

  gamma_vec = (double*) malloc(P->size_of_gamma);
  gsl_matrix_vectorize(gamma_vec, gamma);
  gsl_matrix_free(gamma);
  
  /*PRINTF("cholgam: PDW = %d, LDWORK = %d, 3K = %d\n, s-1 = %d", 1 + (P->s_minus_1 + 1)* P->k_times_d * P->k_times_d, ldwork, 3 * P->k_times_d, P->s_minus_1 );*/

  /* Cholesky factorization of Gamma */
  dwork  = (double*) malloc((size_t)ldwork * sizeof(double));
  mb02gd_("R", "N",  
	  &P->k_times_d, 	/* block size */
	  &P->m_div_k,		/* block_dim(GAMMA) */
	  &P->s_minus_1, 	/* non-zero block superdiagonals */
	  &zero, 		/* previous computations */
	  &P->m_div_k, 	        /* to-be-computed */
	  gamma_vec, /* non-zero part of first block row of GAMMA */
	  &P->k_times_d,      	/* row_dim(gamma) */
	  rb, 			/* packed Cholesky factor */
	  &P->k_times_d_times_s, /* row_dim(rb) */
	  dwork, &ldwork, &info);
  free(dwork);
  free(gamma_vec);

  /* check for errors of mb02gd */
  if (info) { 
    PRINTF("Error: info = %d", info); /* TO BE COMPLETED */
  }
}


/* 
*  JACOBIAN: Computes the Jacobian 
*
*  x_ext  = kron( I_k, [ x; -I_d ] ),
*  params - used for the dimensions n = rowdim(X), 
*           d = coldim(X), k, and n_plus_d = n + d
*  rb     - Choleski factor in a packed form
*  yr     - solution of the system Gamma*yr = r
*  deriv  - Jacobian
*/

void jacobian( gsl_matrix* x_ext, 
	       stls_opt_data* params, 
	       double* rb,
	       gsl_vector* yr,
	       gsl_matrix* deriv )
{
  int i, j, l, k, info;
  const int zero = 0, one = 1;
  gsl_matrix *tmp, *dx_ext, *dgamma, *st;
  gsl_matrix_view submat; 
  gsl_vector_view subvec, res_vec;
  double *dwork, *gamma_vec, *res;

  /* first term of the Jacobian Gamma^{-1/2} kron(a,I_d) */

  /* deriv = kron(a,I_d) */
  gsl_matrix_set_zero(deriv);
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
      submat = gsl_matrix_submatrix(deriv, i*D, j*D, D, D);
      /* assign a(i,j) on the diagonal of submat */
      for (l = 0; l < D; l++)
	gsl_matrix_set(&submat.matrix, l, l, 
		       gsl_matrix_get(P->a, i, j));
    }
  
  /* res = vec(deriv), FORTRAN column-major convention */
  res = malloc(deriv->size2 * deriv->size1 * sizeof(double));
  gsl_matrix_vectorize(res, deriv);

  /* res =  Gamma^{-1/2} * res  */
  dtbtrs_("U", "T", "N", 
	  &P->m_times_d, 
	  &P->k_times_d_times_s_minus_1, 
	  &deriv->size2, 
	  rb, 
	  &P->k_times_d_times_s, 
	  res, 
	  &P->m_times_d, 
	  &info);
  
  /* convert back deriv = res in the C row-major convention */
  gsl_matrix_vec_inv(deriv, res);
  free(res);

  /* second term (naive implementation) */
  dx_ext = gsl_matrix_alloc(SIZE_W, P->k_times_d);
  dgamma = gsl_matrix_alloc(P->k_times_d, 
			   P->k_times_d_times_s);
  tmp  = gsl_matrix_alloc(P->k_times_d, SIZE_W);
  st   = gsl_matrix_alloc(deriv->size1, deriv->size2);
  res  = malloc( P->m_times_d * sizeof(double));
  res_vec = gsl_vector_view_array(res, P->m_times_d);

  for (i = 0; i < N; i++)
    for (j = 0; j < D; j++) {
	/* construct dx_ext */
	gsl_matrix_set_zero(dx_ext);
	for (l = 0; l < P->k; l++)
	  gsl_matrix_set(dx_ext, l*P->n_plus_d+i, l*D+j, 1.0);
	/* form dgamma = d/d x_ij gamma */
	gsl_matrix_set_zero(dgamma);
	for (k = 0; k < S; k++) {
	  submat = gsl_matrix_submatrix(dgamma, 0, k*P->k_times_d, 
					P->k_times_d, P->k_times_d);
	  /* compute tmp = dx_ext' * w_k */
	  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, 
			 dx_ext, P->w->a[k], 0.0, tmp);
	  /* compute submat = tmp * x_ext */
	  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
			 tmp, x_ext, 0.0, &submat.matrix);
	  /* compute tmp = x_ext' * w_k */
	  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, 
			 x_ext, P->w->a[k], 0.0, tmp);
	  /* compute submat = submat  + tmp * dx_ext */
	  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
			 tmp, dx_ext, 1.0, &submat.matrix);
	}
	/* compute st_ij = DGamma * yr */
	subvec = gsl_matrix_column(st, i*D+j);
	tmv_prod(dgamma, S, yr, P->m_div_k, &subvec.vector);
	gsl_vector_memcpy(&res_vec.vector, &subvec.vector);
	/* solve st_ij = Gamma^{-1/2}st_ij */
	dtbtrs_("U", "T", "N", 
		&P->m_times_d, 
		&P->k_times_d_times_s_minus_1, 
		&one, 
		rb, 
		&P->k_times_d_times_s, 
		res,
		&P->m_times_d, 
		&info);
	gsl_vector_memcpy(&subvec.vector, &res_vec.vector);
      }

  free(res);
  gsl_matrix_free(dx_ext);
  gsl_matrix_free(dgamma);
  gsl_matrix_free(tmp);

  /* deriv = deriv - 0.5 * st */
  gsl_matrix_scale(st, 0.5);
  gsl_matrix_sub(deriv, st);

  gsl_matrix_free(st);
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


/* print_arr: print array */
void print_arr(double* a, int n)
{
  int i;

  PRINTF("\n");
  for (i = 0; i < n; i++)
    PRINTF("%f\n",*(a+i));
  PRINTF("\n");
}
