/* stls.c: implementations of the functions from stls.h */

#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>


#include "slra.h"


/* Debug structure for Levenberg-Marquardt algorithm
   typedef struct
   {
   size_t iter;
   double xnorm;
   double fnorm;
   double delta;
   double par;
   gsl_matrix *r;
   gsl_vector *tau;
   gsl_vector *diag;
   gsl_vector *qtf;
   gsl_vector *newton;
   gsl_vector *gradient;
   gsl_vector *x_trial;
   gsl_vector *f_trial;
   gsl_vector *df;
   gsl_vector *sdiag;
   gsl_vector *rptdx;
   gsl_vector *w;
   gsl_vector *work1;
   gsl_permutation * perm;
   }
   lmder_state_t; */


typedef struct
{
  int iter;
  double step;
  double g0norm;
  double pnorm;
  double delta_f;
  double fp0;                   /* f'(0) for f(x-alpha*p) */
  gsl_vector *x0;
  gsl_vector *g0;
  gsl_vector *p;
  /* work space */
  gsl_vector *dx0;
  gsl_vector *dg0;
  gsl_vector *x_alpha;
  gsl_vector *g_alpha;
  /* wrapper function */
  void * wrap;
  /* minimization parameters */
  double rho;
  double sigma;
  double tau1;
  double tau2;
  double tau3;
  int order;
}
vector_bfgs2_state_t;

/*********************/
/* stls: STLS solver */
/*********************/
/* on exit: zero success, otherwise error code: 
   EITER, GSL_ETOLF, GSL_ETOLX, GSL_ETOLG */

static void numer_grad( gsl_vector* x, slra_opt_data_reshaped *P,
			gsl_vector *grad ) {
  gsl_vector *x2;
  double f, f2, dx = 0.00001;
  int i;

  x2 = gsl_vector_alloc(x->size);  
  f = slra_f_reshaped_(x, P);

  for (i = 0; i < x->size; i++) {
    gsl_vector_memcpy(x2, x);
    gsl_vector_set(x2, i, gsl_vector_get(x2, i) + dx);
    f2 = slra_f_reshaped_(x2, P);    
    gsl_vector_set(grad, i, (f2 - f) / dx);
  }
  
  gsl_vector_free(x2);
}

static void comp_meth( gsl_vector* x, slra_opt_data_reshaped *P ) {
  gsl_vector *grad1 = gsl_vector_alloc(P->n * P->d);
  gsl_vector *grad2 = gsl_vector_alloc(P->n * P->d);
  gsl_vector *grad3 = gsl_vector_alloc(P->n * P->d);
  gsl_matrix *jac = gsl_matrix_alloc(P->m * P->d, P->n * P->d);
  gsl_vector *f1 = gsl_vector_alloc(P->m * P->d);
  double val1 = 0, val2 = 0, norm3;

  PRINTF("Compare two methods for derivative:\n");
  slra_fdf_reshaped(x, P, f1, jac);
  slra_fdf_reshaped_(x, P, &val2, grad2);
  numer_grad(x, P, grad3);
  gsl_blas_ddot(f1, f1, &val1);
  gsl_blas_dgemv(CblasTrans, 2.0, jac, f1, 0.0, grad1);

  PRINTF("Cost diff = %f, pseudo grad, grad_new, numer_grad\n",
	 val1 - val2);

  print_arr(grad1->data, grad1->size);
  print_arr(grad2->data, grad2->size);
  print_arr(grad3->data, grad3->size);

  gsl_vector_sub(grad1, grad2);
  gsl_blas_ddot(grad1, grad1, &norm3);
  PRINTF("Deriv diff norm = %f\n", sqrt(norm3));
  
  gsl_vector_free(grad1);
  gsl_vector_free(grad2);
  gsl_vector_free(grad3);
  gsl_vector_free(f1);
  gsl_matrix_free(jac);
}


static int compute_tls( const gsl_matrix *c, const gsl_matrix *perm, gsl_matrix * x ) {
  gsl_matrix * tempc = gsl_matrix_alloc(c->size1, c->size2);
  gsl_matrix_view c_sub_a, c_sub_b;
  int n = x->size1, d = x->size2;
  int status;
    
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, c, perm, 0.0, tempc);
  c_sub_a = gsl_matrix_submatrix(tempc, 0, 0, c->size1, n);
  c_sub_b = gsl_matrix_submatrix(tempc, 0, n, c->size1, d);
    
  status = tls(&c_sub_a.matrix, &c_sub_b.matrix, x);
    
  gsl_matrix_free(tempc);

  if (status) {
    PRINTF("Initial approximation can not be computed: "
        "MB02MD failed with an error info = %d.\n", status);

    return GSL_EINVAL;
  }

  return GSL_SUCCESS;
}


int slra_gsl_optimize( slraFlexCostFunction *F, opt_and_info *opt, gsl_vector* x_vec, gsl_matrix *v ) {
  const gsl_multifit_fdfsolver_type *Tlm[] =
    {gsl_multifit_fdfsolver_lmder, gsl_multifit_fdfsolver_lmsder};

  const gsl_multimin_fdfminimizer_type *Tqn[] = 
    { gsl_multimin_fdfminimizer_vector_bfgs,
      gsl_multimin_fdfminimizer_vector_bfgs2, 
      gsl_multimin_fdfminimizer_conjugate_fr,
      gsl_multimin_fdfminimizer_conjugate_pr };
               
  const gsl_multimin_fminimizer_type *Tnm[] = 
    { gsl_multimin_fminimizer_nmsimplex, gsl_multimin_fminimizer_nmsimplex2, 
      gsl_multimin_fminimizer_nmsimplex2rand };
  
  int gsl_submethod_max[] = { sizeof(Tlm) / sizeof(Tlm[0]),
			  sizeof(Tqn) / sizeof(Tqn[0]),
			  sizeof(Tnm) / sizeof(Tnm[0]) };  
			  
			  
  int status, status_dx, status_grad, k;
  double g_norm, x_norm;


  /* vectorize x row-wise */
  size_t max_ind, min_ind;
  double max_val, min_val, abs_max_val = 0, abs_min_val;
  
  if (opt->method < 0 || 
      opt->method > sizeof(gsl_submethod_max) / sizeof(gsl_submethod_max[0]) || 
      opt->submethod < 0 || 
      opt->submethod > gsl_submethod_max[opt->method]) {
    throw new slraException("Unknown optimization method.\n");   
  }
  

  

  /* LM */
  gsl_multifit_fdfsolver* solverlm;
  gsl_multifit_function_fdf fdflm = { 
    &(slraFlexCostFunction::slra_f_ls), 
    &(slraFlexCostFunction::slra_df_ls), 
    &(slraFlexCostFunction::slra_fdf_ls), F->getM() * F->getD(), F->getN() * F->getD(), F };
  gsl_vector *g;

  /* QN */
  double stepqn = opt->step; /* ??? */
  /*    gsl_multimin_fdfminimizer_vector_bfgs2;*/
  gsl_multimin_fdfminimizer* solverqn;
/*  gsl_multimin_function_fdf fdfqn = { 
    &slra_f_reshaped_,
    &slra_df_reshaped_, 
    &slra_fdf_reshaped_, P->n * P->d, P }; */

  gsl_multimin_function_fdf fdfqn = { 
    &(slraFlexCostFunction::slra_f),
    &(slraFlexCostFunction::slra_df), 
    &(slraFlexCostFunction::slra_fdf), F->getN() * F->getD(), F };


  /* NM */
  double size;
  gsl_vector *stepnm;
  gsl_multimin_fminimizer* solvernm;
  /*gsl_multimin_function fnm = { &slra_f_reshaped_, P->n * P->d, P };*/
  gsl_multimin_function fnm = { &slra_f_reshaped_, F->getN() * F->getD(), F };
  
  /* comp_meth(&x_vec.vector, &params); */

  /* initialize the optimization method */
  switch (opt->method) {
  case SLRA_OPT_METHOD_LM: /* LM */
    solverlm = gsl_multifit_fdfsolver_alloc(Tlm[opt->submethod], 
					    F->getM() * F->getD(), F->getN() * F->getD());
    gsl_multifit_fdfsolver_set(solverlm, &fdflm, x_vec);
    g = gsl_vector_alloc(F->getN() * F->getD());
    break;
  case SLRA_OPT_METHOD_QN: /* QN */
    solverqn = gsl_multimin_fdfminimizer_alloc( Tqn[opt->submethod], 
						F->getN() * F->getD() );
    gsl_multimin_fdfminimizer_set(solverqn, &fdfqn, x_vec, 
				  stepqn, opt->tol); 
    status_dx = GSL_CONTINUE;  
    break;
  case SLRA_OPT_METHOD_NM: /* NM */
    solvernm = gsl_multimin_fminimizer_alloc( Tnm[opt->submethod], F->getN() * F->getD() );
    stepnm = gsl_vector_alloc( F->getN() * F->getD() );
    gsl_vector_set_all(stepnm, opt->step); 
    gsl_multimin_fminimizer_set( solvernm, &fnm, x_vec, stepnm );
    break;
  }

  /* optimization loop */
  if (opt->disp == SLRA_OPT_DISP_FINAL || opt->disp == SLRA_OPT_DISP_ITER) {
    PRINTF("STLS optimization:\n");
  }
    
  status = GSL_SUCCESS;  
  status_dx = GSL_CONTINUE;
  status_grad = GSL_CONTINUE;  
  opt->iter = 0;
  
  
  if (opt->method == SLRA_OPT_METHOD_LM && opt->disp == SLRA_OPT_DISP_ITER) {
    gsl_blas_ddot(solverlm->f, solverlm->f, &opt->fmin);

    gsl_multifit_gradient(solverlm->J, solverlm->f, g);	
    x_norm = gsl_blas_dnrm2(solverlm->x);
    g_norm = gsl_blas_dnrm2(g);
    PRINTF("  0: f0 = %16.8f,  ||f0'|| = %16.8f,  ||x|| = %10.8f\n", opt->fmin, g_norm, x_norm);
  }
  
  
  if (opt->submethod == SLRA_OPT_SUBMETHOD_QN_BFGS2 && opt->disp == SLRA_OPT_DISP_ITER) {
     vector_bfgs2_state_t *st = (vector_bfgs2_state_t *) solverqn->state;
      
     PRINTF("State info: pnorm = %16.8f,  g0norm = %16.8f, fp0 = %16.8f\n",
         st->pnorm, st->g0norm, st->fp0);
          
          
  }
  
  while (status_dx == GSL_CONTINUE && 
	 status_grad == GSL_CONTINUE &&
	 status == GSL_SUCCESS &&
	 opt->iter < opt->maxiter) {
    /* print_vec(solverlm->x); */
    opt->iter++;
    switch (opt->method) {
    case SLRA_OPT_METHOD_LM: /* Levenberge-Marquardt */
      status = gsl_multifit_fdfsolver_iterate( solverlm );
      /* check for convergence problems */
      if (status == GSL_ETOLF || 
	  status == GSL_ETOLX || 
	  status == GSL_ETOLG) {
	break; /* <- THIS IS WRONG */
      }
      /* check the convergence criteria */
      status_dx = gsl_multifit_test_delta(solverlm->dx, 
					  solverlm->x, 
					  opt->epsabs, 
					  opt->epsrel);
      gsl_multifit_gradient(solverlm->J, solverlm->f, g);
      status_grad = gsl_multifit_test_gradient(g, opt->epsgrad);
      /* print information */
      if (opt->disp == SLRA_OPT_DISP_ITER) {
	gsl_blas_ddot(solverlm->f, solverlm->f, &opt->fmin);
	
	x_norm = gsl_blas_dnrm2(solverlm->x);
	g_norm = gsl_blas_dnrm2(g);
	
	PRINTF("%3u: f0 = %16.8f,  ||f0'|| = %16.8f,  ||x|| = %10.8f\n",
	       opt->iter, opt->fmin, g_norm, x_norm);
      }
      break;
    case SLRA_OPT_METHOD_QN:
      status = gsl_multimin_fdfminimizer_iterate( solverqn );
      if (status == GSL_ENOPROG) {
        if (opt->submethod == SLRA_OPT_SUBMETHOD_QN_BFGS2 && opt->disp == SLRA_OPT_DISP_ITER) {
          vector_bfgs2_state_t *st = (vector_bfgs2_state_t *)solverqn->state;
          
          PRINTF("No progress: pnorm = %16.8f,  g0norm = %16.8f, fp0 = %16.8f\n",
              st->pnorm, st->g0norm, st->fp0);
          
          
        }
	break; /* <- THIS IS WRONG */
      }

      /* check the convergence criteria */
      status_grad = gsl_multimin_test_gradient(
		    gsl_multimin_fdfminimizer_gradient(solverqn), 
		    opt->epsgrad );
		    
     status_dx = gsl_multifit_test_delta(solverqn->dx, 
                                         solverqn->x, 
					 opt->epsabs, 
					 opt->epsrel);  		    
      if (opt->disp == SLRA_OPT_DISP_ITER) {
	opt->fmin = gsl_multimin_fdfminimizer_minimum( solverqn );
	x_norm = gsl_blas_dnrm2(solverqn->x);
	g_norm = gsl_blas_dnrm2(solverqn->gradient);
	PRINTF("%3u: f0 = %16.8f,  ||f0'|| = %16.8f,  ||x|| = %10.8f\n", 
	       opt->iter, opt->fmin, g_norm, x_norm);
      }
      break;
    case SLRA_OPT_METHOD_NM:
      status = gsl_multimin_fminimizer_iterate( solvernm );
      /* check the convergence criteria */
      size = gsl_multimin_fminimizer_size( solvernm );
      status_dx = gsl_multimin_test_size( size, opt->epsx );
      /* print information */
      if (opt->disp == SLRA_OPT_DISP_ITER) {
	opt->fmin = gsl_multimin_fminimizer_minimum( solvernm );
	x_norm = gsl_blas_dnrm2(solvernm->x);

	PRINTF("%3u: f0 = %16.8f,  ||x|| = %10.8f\n", 
	       opt->iter, opt->fmin, g_norm, x_norm);
      }
      break;
    }
  } 
  if (opt->iter >= opt->maxiter) {
    status = EITER;
  }

  switch (opt->method) {
  case  SLRA_OPT_METHOD_LM:
    /* return the results */
    gsl_vector_memcpy(x_vec, solverlm->x);
    gsl_multifit_covar(solverlm->J, opt->epsrel, v); /* ??? Different eps */
    /* assign the opt output fields */
    gsl_blas_ddot(solverlm->f, solverlm->f, &opt->fmin);
    break;
  case SLRA_OPT_METHOD_QN:
    gsl_vector_memcpy(x_vec, solverqn->x);

    opt->fmin = solverqn->f;
    break;
  case SLRA_OPT_METHOD_NM:
    /* return the results */
    gsl_vector_memcpy(x_vec, solvernm->x);
    /* gsl_multifit_covar( J??, opt->epsrel, v); */
    /* assign the opt output fields */
    opt->fmin = solvernm->fval;
    break;
  }
  
/*  if (opt->disp == SLRA_OPT_DISP_ITER) {
     opt->chol_time =  ((double)P->chol_time / CLOCKS_PER_SEC) / P->chol_count;
  }*/

  
  /* print exit information */  
  if (opt->disp != SLRA_OPT_DISP_OFF) { /* unless "off" */
    switch (status) {
    case EITER: 
      PRINTF("STLS optimization terminated by reaching the maximum number " 
	     "of iterations.\nThe result could be far from optimal.\n");
      break;
    case GSL_ETOLF:
      PRINTF("Lack of convergence: progress in function value < machine EPS.\n");
      break;
    case GSL_ETOLX:
      PRINTF("Lack of convergence: change in parameters < machine EPS.\n");
      break;
    case GSL_ETOLG:
      PRINTF("Lack of convergence: change in gradient < machine EPS.\n");
      break;
    case GSL_ENOPROG:
      PRINTF("Possible lack of convergence: no progress.\n");
      break;
    }
    if (!status && (opt->disp == SLRA_OPT_DISP_FINAL || 
                    opt->disp == SLRA_OPT_DISP_ITER)) { 
      if (status_grad == GSL_CONTINUE) {
	PRINTF("Optimization terminated by reaching the convergence " 
	       "tolerance for X.\n");
      } else if (status_dx == GSL_CONTINUE) {
	PRINTF("Optimization terminated by reaching the convergence " 
	       "tolerance for the gradient.\n");
      } else {
	PRINTF("Optimization terminated by reaching the convergence " 
	       "tolerance for both X and the gradient.\n"); 
      }
    }
  }

  /* Cleanup  */
  switch (opt->method) {
  case SLRA_OPT_METHOD_LM: /* LM */
    gsl_multifit_fdfsolver_free(solverlm);
    gsl_vector_free(g);
    break;
  case SLRA_OPT_METHOD_QN: /* QN */
    gsl_multimin_fdfminimizer_free(solverqn);
    break;
  case SLRA_OPT_METHOD_NM: /* NM */
    gsl_multimin_fminimizer_free(solvernm);
    gsl_vector_free(stepnm);
    break;
  }


  return GSL_SUCCESS; /* <- correct with status */
}


static void checkAndComputeX( slraFlexCostFunction * F, gsl_matrix *x, opt_and_info* opt, int x_given ) {
  if (x->size1 + x->size2 != F->getNplusD()) {
    throw new slraException("Initial approximation doesn't conform to the structure" 
                            "specification.\n");   
  }      

  if (!x_given) {  /* compute default initial approximation */
    if (opt->disp == SLRA_OPT_DISP_ITER) {
      PRINTF("X not given, computing TLS initial approximation.\n");
    }
    if (compute_tls(F->getSMatr(), F->getPerm(), x) ==  GSL_EINVAL) {
      throw new slraException("Error while computing initial approximation.\n");   
    }
  }
}


int slra(gsl_vector* p, data_struct* s, int r, gsl_matrix* x,
         gsl_matrix* v, opt_and_info* opt, int x_given, int compute_ph,
         gsl_matrix *perm ) {
  slraFlexCostFunction * myCostFun = NULL;
  int res = GSL_SUCCESS;

  try { 
    myCostFun =  new slraFlexCostFunction(slraFlexStructure(s, p->size), r, p, opt, perm);
    checkAndComputeX(myCostFun, x, opt, x_given);
  
    time_t t_b = clock();
    gsl_vector_view x_vec = gsl_vector_view_array(x->data, x->size1 * x->size2);
    int status = slra_gsl_optimize(myCostFun, opt, &(x_vec.vector), v);
    opt->time = (double) (clock() - t_b) / (double) CLOCKS_PER_SEC;

    if (compute_ph) {
      myCostFun->computeCorrection(p, &(x_vec.vector));
    }
  } catch ( slraException *e ) {
    res = GSL_EINVAL;
    PRINTF(e->getMessage());
    delete e;
  }

  if (myCostFun != NULL) {
    delete myCostFun;
  }

  return res;
}

/*#define P ((slra_opt_data*) params)*/

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
  if ( m < n + l ) {
    ldwork = GSL_MAX(m * (n + l) + GSL_MAX(3 * m + n + l, 5 * m), 3 * l);
  } else {
    ldwork = GSL_MAX(3*(n+l)+m, 5*(n+l));
  }
  
  /* c = [a b] in FORTRAN format */
  c = (double *)malloc( ldc * (n+l) * sizeof(double));
  gsl_matrix_vectorize(c, a);
  gsl_matrix_vectorize(c + m*n, b);

  xa = (double *)malloc( n * l * sizeof(double));
  s  = (double *)malloc( (n+l) * sizeof(double));
  dwork = (double *)malloc( ldwork * sizeof(double));
  iwork = (int *)malloc( l * sizeof(int));

  /* Call MB02MD */
  mb02md_("B", &m, &n, &l, &r, c, &ldc, s, xa, &n, 
	  &tol, iwork, dwork, &ldwork, &iwarn, &info);

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
