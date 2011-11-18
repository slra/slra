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




int stls(gsl_matrix* a, gsl_matrix* b, data_struct* s, 
         gsl_matrix* x, gsl_matrix* v, opt_and_info* opt, gsl_vector *p) {
  int status;
  char method = 'l';
  int status_dx, status_grad, k;
  double g_norm, x_norm;
  int m,n,d;
  int np, n_plus_d;
  time_t t_b;
  //stls_opt_data_old params;
  stls_opt_data_reshaped params;


  gsl_matrix *c;
  gsl_matrix_view c_sub_a, c_sub_b;
  

  t_b = clock();
  /* constants */
  n = x->size1;
  d = x->size2;
  
  if (p != NULL) {
    if (a != NULL && s->m != a->size1) {
      PRINTF("Error in structure specification: m is not correct in s.\n");   
      return GSL_EINVAL;
    }
    
    if (check_and_adjust_parameters(s, &np, &n_plus_d) != GSL_SUCCESS) {
      PRINTF("Error in structure specification: incorrect number of columns/rows in a block.\n");   
      return GSL_EINVAL;
    }
    
    if (n+d != n_plus_d) {
      PRINTF("Error in structure specification: total number of columns.\n");   
      return GSL_EINVAL;
    }

    
    if (np != p->size) {
      PRINTF("Error in structure specification: incorrect length of the parameter vector.\n");
      return GSL_EINVAL;
    }
  
    c = gsl_matrix_alloc(s->m, n_plus_d);
    stls_fill_matrix_from_p(c, s, p);
    c_sub_a = gsl_matrix_submatrix(c, 0, 0, s->m, n);
    c_sub_b = gsl_matrix_submatrix(c, 0, n, s->m, d);
  } else {
    if (a == NULL) {
      PRINTF("Both a and p are NULL specified \n");
      return GSL_EINVAL;
    }
  }


  if (a == NULL) {
    m = s->m;
    a = &c_sub_a.matrix;
    b = &c_sub_b.matrix;
  } else {
    m = a->size1;
    
    gsl_matrix_sub(&c_sub_a.matrix, a);
    gsl_matrix_sub(&c_sub_b.matrix, b);
    
    PRINTF("a and b don't coincide with a and b computed from p\n");
    double minc, maxc;
    
    gsl_matrix_minmax(c, &minc,&maxc);
    
    if (minc != 0.0 || maxc != 0.0) {
      PRINTF("a and b don't coincide with a and b computed from p\n");
    }
    
    /* !! Comparison of old parameters and new */
    
  } 


  //allocate_and_prepare_data_old(a, b, s, opt,  &params);
  allocate_and_prepare_data_reshaped(a, b, s, opt, &params);

  /* LM */
  const gsl_multifit_fdfsolver_type *Tlm
    = gsl_multifit_fdfsolver_lmder;
  gsl_multifit_fdfsolver* solverlm;
//  gsl_multifit_function_fdf fdflm = { &stls_f, &stls_df, &stls_fdf, m * d, n * d, &params};
  gsl_multifit_function_fdf fdflm = { &stls_f_reshaped, &stls_df_reshaped, &stls_fdf_reshaped, m * d, n * d, &params};
  gsl_vector *g;

  /* QN */
  double stepqn = 0.001; /* ??? */
  const gsl_multimin_fdfminimizer_type *Tqn
    = gsl_multimin_fdfminimizer_vector_bfgs;
  gsl_multimin_fdfminimizer* solverqn;
  gsl_multimin_function_fdf fdfqn;
//  gsl_multimin_function_fdf fdfqn = { &stls_f_, &stls_df_, &stls_fdf_, P->n_times_d, P };

  /* NM */
  double size;
  gsl_vector *stepnm;
  const gsl_multimin_fminimizer_type *Tnm
    = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer* solvernm;
  gsl_multimin_function fnm = { &stls_f_reshaped_, n * d, &params };

  /* vectorize x row-wise */
  gsl_vector_view x_vec;
  x_vec = gsl_vector_view_array(x->data, n * d);

  /* initialize the optimization method */
  switch (method) {
  case 'l': /* LM */
    solverlm = gsl_multifit_fdfsolver_alloc(Tlm, m * d, n * d);
    gsl_multifit_fdfsolver_set(solverlm, &fdflm, &x_vec.vector);
    g = gsl_vector_alloc(n * d);
    break;
  case 'q': /* QN * /

    solverqn = gsl_multimin_fdfminimizer_alloc(Tqn, P->n_times_d);
    gsl_multimin_fdfminimizer_set(solverqn, &fdfqn, &x_vec.vector, stepqn, opt->epsabs);
    status_dx = GSL_CONTINUE;  */
    break;
  case 'n': /* NM */
    solvernm = gsl_multimin_fminimizer_alloc( Tnm, n * d );
    stepnm = gsl_vector_alloc( n * d );
    gsl_vector_set_all( stepnm, 0.001 ); /* ??? */
    gsl_multimin_fminimizer_set( solvernm, &fnm, &x_vec.vector, stepnm );
    break;
  default:
    ; /* error */
  }

  /* optimization loop */
  if (opt->disp == 2 || opt->disp == 3)
    PRINTF("STLS optimization:\n");
    
    
  status_dx = GSL_CONTINUE;
  status_grad = GSL_CONTINUE;  
  opt->iter = 0;
  
  
  while (status_dx == GSL_CONTINUE && status_grad == GSL_CONTINUE && opt->iter < opt->maxiter) {
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
	
	x_norm = gsl_blas_dnrm2(solverlm->x);
	g_norm = gsl_blas_dnrm2(g);
	PRINTF("%3u: f0 = %16.8f,  ||f0'|| = %16.8f,  ||x|| = %13.5f\n", 
	       opt->iter, opt->fmin, g_norm, x_norm);
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
  } 
  if (opt->iter >= opt->maxiter) {
    status = EITER;
  }
  opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;

  switch (method) {
  case  'l':
    /* return the results */
    gsl_vector_memcpy(&x_vec.vector, solverlm->x);
    gsl_multifit_covar(solverlm->J, opt->epsrel, v);
    /* assign the opt output fields */
    gsl_blas_ddot(solverlm->f, solverlm->f, &opt->fmin);
    break;
  case 'q':
    break;
  case 'n':
    /* return the results */
    gsl_vector_memcpy(&x_vec.vector, solvernm->x);
    /* gsl_multifit_covar( J??, opt->epsrel, v); */
    /* assign the opt output fields */
    opt->fmin = solvernm->fval;
    break;
  }


  if (p != NULL) {
    if (opt->corr) {
      stls_correction_reshaped(p, s, &params, &(x_vec.vector));
    }
    gsl_matrix_free(c);
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

  return GSL_SUCCESS; /* <- correct with status */
}



/*#define P ((stls_opt_data*) params)*/





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


