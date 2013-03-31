#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "slra.h"

OptimizationOptions::OptimizationOptions() : disp(SLRA_DEF_disp), 
    method(SLRA_DEF_method),
    submethod(SLRA_DEF_submethod),  maxiter(SLRA_DEF_maxiter),
    epsabs(SLRA_DEF_epsabs), epsrel(SLRA_DEF_epsrel), 
    epsgrad(SLRA_DEF_epsgrad), epsx(SLRA_DEF_epsx), maxx(SLRA_DEF_maxx),
    step(SLRA_DEF_step), tol(SLRA_DEF_tol), reggamma(SLRA_DEF_reggamma),
    ls_correction(SLRA_DEF_ls_correction) {
}

int OptimizationOptions::gslOptimize( OptFunction *F, 
                                      gsl_vector* x_vec, gsl_matrix *v ) {
  const gsl_multifit_fdfsolver_type *Tlm[] =
    { gsl_multifit_fdfsolver_lmder, gsl_multifit_fdfsolver_lmsder };
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
  
  if (this->method < 0 || 
      this->method > sizeof(gsl_submethod_max)/sizeof(gsl_submethod_max[0]) || 
      this->submethod < 0 || 
      this->submethod > gsl_submethod_max[this->method]) {
    throw new Exception("Unknown optimization method.\n");   
  }
  
  if (this->maxiter < 0 || this->maxiter > 5000) {
    throw new Exception("opt.maxiter should be in [0;5000].\n");   
  }

  /* LM */
  gsl_multifit_fdfsolver* solverlm;
  gsl_multifit_function_fdf fdflm = { &(F->_f_ls),  &(F->_df_ls), &(F->_fdf_ls), 
                                       F->getNsq(), F->getNvar(), F };
  gsl_vector *g;

  /* QN */
  double stepqn = this->step; 
  gsl_multimin_fdfminimizer* solverqn;
  gsl_multimin_function_fdf fdfqn = { 
    &(F->_f), &(F->_df), &(F->_fdf), F->getNvar(), F };

  /* NM */
  double size;
  gsl_vector *stepnm;
  gsl_multimin_fminimizer* solvernm;
  gsl_multimin_function fnm = { &(F->_f), F->getNvar(), F };

  /* initialize the optimization method */
  switch (this->method) {
  case SLRA_OPT_METHOD_LM: /* LM */
    solverlm = gsl_multifit_fdfsolver_alloc(Tlm[this->submethod], 
                   F->getNsq(), F->getNvar());
    gsl_multifit_fdfsolver_set(solverlm, &fdflm, x_vec);
    g = gsl_vector_alloc(F->getNvar());
    break;
  case SLRA_OPT_METHOD_QN: /* QN */
    solverqn = gsl_multimin_fdfminimizer_alloc(Tqn[this->submethod], 
						F->getNvar() );
    gsl_multimin_fdfminimizer_set(solverqn, &fdfqn, x_vec, 
				  stepqn, this->tol); 
    status_dx = GSL_CONTINUE;  
    break;
  case SLRA_OPT_METHOD_NM: /* NM */
    solvernm = gsl_multimin_fminimizer_alloc(Tnm[this->submethod], F->getNvar());
    stepnm = gsl_vector_alloc(F->getNvar());
    gsl_vector_set_all(stepnm, this->step); 
    gsl_multimin_fminimizer_set( solvernm, &fnm, x_vec, stepnm );
    break;
  }

  /* optimization loop */
  Log::lprintf(Log::LOG_LEVEL_FINAL, "SLRA optimization:\n");
    
  status = GSL_SUCCESS;  
  status_dx = GSL_CONTINUE;
  status_grad = GSL_CONTINUE;  
  this->iter = 0;
  
  if ((this->method == SLRA_OPT_METHOD_LM || 
       this->method == SLRA_OPT_METHOD_QN) && 
      Log::getMaxLevel() >= Log::LOG_LEVEL_ITER) {
    if (this->method == SLRA_OPT_METHOD_LM) { 
      gsl_blas_ddot(solverlm->f, solverlm->f, &this->fmin);
      gsl_multifit_gradient(solverlm->J, solverlm->f, g);	
      x_norm = gsl_blas_dnrm2(solverlm->x);
      g_norm = gsl_blas_dnrm2(g);
    } else {
      this->fmin = gsl_multimin_fdfminimizer_minimum(solverqn);
      x_norm = gsl_blas_dnrm2(solverqn->x);
      g_norm = gsl_blas_dnrm2(solverqn->gradient);
    }
    Log::lprintf("  0: f0 = %15.10e,  ||f0'|| = %15.7e,  ||x|| = %10.8f\n",
                this->fmin, sqrt(2) * g_norm, x_norm);
  }

  while (status_dx == GSL_CONTINUE && 
	 status_grad == GSL_CONTINUE &&
	 status == GSL_SUCCESS &&
	 this->iter < this->maxiter) {
  	if (this->method == SLRA_OPT_METHOD_LM && this->maxx > 0) {
  	  if (gsl_vector_max(solverlm->x) > this->maxx || gsl_vector_min(solverlm->x) < -this->maxx){
  	    break;
	    }
	  }

    /* print_vec(solverlm->x); */
    this->iter++;
    switch (this->method) {
    case SLRA_OPT_METHOD_LM: /* Levenberge-Marquardt */

      status = gsl_multifit_fdfsolver_iterate(solverlm);
      /* check for convergence problems */
      if (status == GSL_ETOLF || 
          status == GSL_ETOLX || 
          status == GSL_ETOLG) {
        break; /* <- THIS IS WRONG */
      }
      /* check the convergence criteria */
      if (this->epsabs != 0 || this->epsrel != 0) {
        status_dx = gsl_multifit_test_delta(solverlm->dx, solverlm->x, 
	  				  this->epsabs, this->epsrel);
	  	} else {
	  	  status_dx = GSL_CONTINUE;
	  	}
      gsl_multifit_gradient(solverlm->J, solverlm->f, g);
      status_grad = gsl_multifit_test_gradient(g, this->epsgrad);
      /* print information */
      if (Log::getMaxLevel() >= Log::LOG_LEVEL_ITER) {
        gsl_blas_ddot(solverlm->f, solverlm->f, &this->fmin);
        x_norm = gsl_blas_dnrm2(solverlm->x);
        g_norm = gsl_blas_dnrm2(g);
        Log::lprintf("%3u: f0 = %15.10e,  ||f0'|| = %15.7e,  ||x|| = %10.8f\n",
                     this->iter, this->fmin, sqrt(2) * g_norm, x_norm);
      }  
      break;
    case SLRA_OPT_METHOD_QN:
      status = gsl_multimin_fdfminimizer_iterate( solverqn );
/*      if (status == GSL_ENOPROG) {
	break; /* <- THIS IS WRONG * /
      }*/

      /* check the convergence criteria */
      status_grad = gsl_multimin_test_gradient(
          gsl_multimin_fdfminimizer_gradient(solverqn), this->epsgrad );
		    
      status_dx = gsl_multifit_test_delta(solverqn->dx, solverqn->x, 
	 				 this->epsabs, this->epsrel);  		    
      if (Log::getMaxLevel() >= Log::LOG_LEVEL_ITER) {
        this->fmin = gsl_multimin_fdfminimizer_minimum(solverqn);
        x_norm = gsl_blas_dnrm2(solverqn->x);
        g_norm = gsl_blas_dnrm2(solverqn->gradient);
        Log::lprintf("%3u: f0 = %16.11f,  ||f0'|| = %16.8f,  ||x|| = %10.8f\n", 
                    this->iter, this->fmin, g_norm, x_norm);
      }
      break;
    case SLRA_OPT_METHOD_NM:
      status = gsl_multimin_fminimizer_iterate( solvernm );
      /* check the convergence criteria */
      size = gsl_multimin_fminimizer_size( solvernm );
      status_dx = gsl_multimin_test_size( size, this->epsx );
      /* print information */
      if (Log::getMaxLevel() >= Log::LOG_LEVEL_ITER) {
        this->fmin = gsl_multimin_fminimizer_minimum( solvernm );
        x_norm = gsl_blas_dnrm2(solvernm->x);
        Log::lprintf("%3u: f0 = %15.10e,  ||x|| = %9.7e\n", 
	            this->iter, this->fmin, g_norm, x_norm);
      }
      break;
    }
  } 
  if (this->iter >= this->maxiter) {
    status = EITER;
  }

  switch (this->method) {
  case  SLRA_OPT_METHOD_LM:
    /* return the results */
    gsl_vector_memcpy(x_vec, solverlm->x);
    if (v != NULL) {
      gsl_multifit_covar(solverlm->J, this->epscov, v); /* ??? Different eps */
    }
    /* assign the opt output fields */
    gsl_blas_ddot(solverlm->f, solverlm->f, &this->fmin);
    break;
  case SLRA_OPT_METHOD_QN:
    gsl_vector_memcpy(x_vec, solverqn->x);
    this->fmin = solverqn->f;
    break;
  case SLRA_OPT_METHOD_NM:
    /* return the results */
    gsl_vector_memcpy(x_vec, solvernm->x);
    /* gsl_multifit_covar( J??, this->epsrel, v); */
    /* assign the opt output fields */
    this->fmin = solvernm->fval;
    break;
  }
  
  /* print exit information */  
  if (Log::getMaxLevel() >= Log::LOG_LEVEL_FINAL) { /* unless "off" */
    switch (status) {
    case EITER: 
      Log::lprintf("SLRA optimization terminated by reaching " 
                  "the maximum number of iterations.\n" 
                  "The result could be far from optimal.\n");
      break;
    case GSL_ETOLF:
      Log::lprintf("Lack of convergence: "
                  "progress in function value < machine EPS.\n");
      break;
    case GSL_ETOLX:
      Log::lprintf("Lack of convergence: "
                  "change in parameters < machine EPS.\n");
      break;
    case GSL_ETOLG:
      Log::lprintf("Lack of convergence: "
                  "change in gradient < machine EPS.\n");
      break;
    case GSL_ENOPROG:
      Log::lprintf("Possible lack of convergence: no progress.\n");
      break;
    }
    
    if (status_grad != GSL_CONTINUE && status_dx != GSL_CONTINUE) {
      Log::lprintf("Optimization terminated by reaching the convergence "
                  "tolerance for both X and the gradient.\n"); 
    
    } else {
      if (status_grad != GSL_CONTINUE) {
	Log::lprintf("Optimization terminated by reaching the convergence "
	            "tolerance for the gradient.\n");
      } else {
        Log::lprintf("Optimization terminated by reaching the convergence "
                    "tolerance for X.\n");
      }
    }
  }

  /* Cleanup  */
  switch (this->method) {
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



