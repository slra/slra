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
  int status;
  int m,n,d;
  char method = 'l';
  int status_dx, status_grad, k;
  double g_norm;
  //stls_opt_data_old params;
  stls_opt_data_reshaped params;
  

  
  
/*  printf("P = %p, sizeof(stls_opt_data) = %i\n", P, sizeof(stls_opt_data));*/


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
  
  
  //allocate_and_prepare_data_old(a, b, s, P);

  allocate_and_prepare_data_reshaped(a, b, s, opt, &params);
  //allocate_and_prepare_data_old(a, b, s, opt,  &params);


  /* vectorize x row-wise */
  x_vec = gsl_vector_view_array(x->data, n * d);

  /* initiaalize the optimization method */
  switch (method) {
  case 'l': /* LM */
  /** /
    fdflm.f      = &stls_f;
    fdflm.df     = &stls_df;
    fdflm.fdf    = &stls_fdf; /**/

    fdflm.f      = &stls_f_reshaped;
    fdflm.df     = &stls_df_reshaped;
    fdflm.fdf    = &stls_fdf_reshaped; /**/
    
    fdflm.n      = m * d;
    fdflm.p      = n * d;
    fdflm.params = &params;
    
    solverlm = gsl_multifit_fdfsolver_alloc(Tlm, m * d, n * d);
    gsl_multifit_fdfsolver_set(solverlm, &fdflm, &x_vec.vector);
    g = gsl_vector_alloc(n * d);
    break;
  case 'q': /* QN 
    fdfqn.f      = &stls_f_;
    fdfqn.df     = &stls_df_;
    fdfqn.fdf    = &stls_fdf_;
    fdfqn.n      = P->n_times_d;
    fdfqn.params = P;

    solverqn = gsl_multimin_fdfminimizer_alloc(Tqn, P->n_times_d);
    gsl_multimin_fdfminimizer_set(solverqn, &fdfqn, &x_vec.vector, stepqn, opt->epsabs);
    status_dx = GSL_CONTINUE; */
    break;
  case 'n': /* NM */
    fnm.f = &stls_f_reshaped_;
    fnm.n = n * d;
    fnm.params = &params;

    solvernm = gsl_multimin_fminimizer_alloc( Tnm, n * d );
    stepnm = gsl_vector_alloc( n * d );
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

   free_memory_reshaped(&params);

// free_memory_old(&params);





/*  free(P);*/

  return GSL_SUCCESS; /* <- correct with status */
}



/*#define P ((stls_opt_data*) params)*/



/*
* tmv_prod: block-Toeplitz matrix T times vector v
*
* tt - storage for [t_s-1' ... t_1' t_0 t_1 ... t_s-1].
* s - number of blocks in t, t = [t_0 ... t_s-1]
* s_1 = s - 1;  m = (int) v->size1 / tt->size1
* p - result
*/ 
void tmv_prod_new( gsl_matrix *tt, int s, gsl_vector* v, int m, 
	      gsl_vector* p)
{
  int i, imax, temp, s_1 = s - 1;
  int row_lim = GSL_MIN(s_1, m/2);
  gsl_vector_view subv, subp; 	/* subvector of v and p */

  int TM = tt->size1; 		/* = block size */

  gsl_matrix_view submat, source;


  /* form tt */


/*  PRINTF("s = %d, m = %d, row_lim = %d", s, m, row_lim);*/
  /* construct p = T*v */
  gsl_vector_set_zero(p);

  /* beginning and end parts of the product p */
  for (i = 0; i < row_lim; i++) {
    temp = GSL_MIN(s+i, m)*TM;
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
    submat = gsl_matrix_submatrix(tt, 0, (s+i)*TM -temp, TM, temp);    
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







/*
* tmv_prod: block-Toeplitz matrix T times vector v
*
* t - nonzero part of the first block row of T
* s - number of blocks in t, t = [t_0 ... t_s-1]
* s_1 = s - 1;  m = (int) v->size1 / t->size1
* p - result
* 
*/ 

void tmv_prod(gsl_matrix* t, int s, gsl_vector* v, int m, 
	      gsl_vector* p)
{
  gsl_matrix *tt;
  int i, imax, temp, s_1 = s - 1;
  gsl_vector_view subv, subp; 	/* subvector of v and p */

  int TM = t->size1; 		/* = block size */
  int TN = (t->size2);		/* = s(block size) */


  /* tt - storage for [t_s-1' ... t_1' t_0 t_1 ... t_s-1]. Should be t->size1* (2 * t->size2 - t->size1). */

  gsl_matrix_view submat, source;
  


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

  tmv_prod_new(tt, s,v, m, p);

  free(tt);
  
}




 /* find n_d = n+d = sum_{l=1}^q n_l and w->s */
int get_bandwidth_from_structure( const data_struct* s ) {
  int l, max_nl = 1;

  for (l = 0; l < s->q; l++) {
    if ((s->a[l].type == 'T' || s->a[l].type == 'H')) {
      max_nl = mymax(s->a[l].ncol / s->a[l].nb, max_nl);
    }
  }
  
  return max_nl;
}




/* s2w: finds the covariance matrices w from the data structure */
int s2w(const data_struct* s, w_data* w, int n_d, int blocked )
{
  int k, l, i, offset, imax, sum_nl;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;
  char err_msg[70];
  int rep;
  int size_wk;
  
 
  w->s = get_bandwidth_from_structure(s);
  
  if (blocked) {
    rep = s->k;
    size_wk = s->k * n_d;
  } else {
    rep = 1;
    size_wk = n_d;
  }

  w->a = (gsl_matrix**) malloc(w->s * sizeof(gsl_matrix *));
  zk   = gsl_matrix_alloc(n_d, n_d);
  /* construct w */
  for (k = 0; k < w->s; k++) { 
    gsl_matrix_set_zero(zk);
    for (l = sum_nl = 0; l < s->q; sum_nl += s->a[l++].ncol) { 
      zkl = gsl_matrix_submatrix(zk, sum_nl, sum_nl, s->a[l].ncol, s->a[l].ncol); 
      switch (s->a[l].type) {
      case 'T': case 'H':
	offset = s->a[l].nb * k;
	imax   = s->a[l].ncol  - offset;
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
    w->a[k] = gsl_matrix_calloc(size_wk, size_wk);
    /* w->a[k] = kron(Ik, zk) */
    for (i = 0; i < rep; i++) {
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
  c = malloc( ldc * (n+l) * sizeof(double));
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

  for (j = 0; j < m->size2; j++)
    for (i = 0; i < m->size1; i++)
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
