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
    ls_correction(SLRA_DEF_ls_correction), avoid_xi(SLRA_DEF_avoid_xi) {
}

void OptimizationOptions::str2Method( const char *str )  {
  char meth_codes[] = "lqnp", 
       sm_codes_lm[] = "ls", sm_codes_qn[] = "b2pf", sm_codes_nm[] = "n2r", sm_codes_lmpinv[] = "su";
  char *submeth_codes[] = { sm_codes_lm, sm_codes_qn, sm_codes_nm, sm_codes_lmpinv };

  size_t submeth_codes_max[] = { 
    sizeof(sm_codes_lm) / sizeof(sm_codes_lm[0]) - 1, 
    sizeof(sm_codes_qn) / sizeof(sm_codes_qn[0]) - 1, 
    sizeof(sm_codes_nm) / sizeof(sm_codes_nm[0]) - 1,
    sizeof(sm_codes_lmpinv) / sizeof(sm_codes_lmpinv[0]) - 1
  };
  size_t meth_code_max = sizeof(submeth_codes_max) / sizeof(submeth_codes_max[0]);
  long i;
  
  if (str[0] == 0) {
    return;
  }
  for (i = meth_code_max - 1; i >= 0; --i)  {
    if (str[0] == meth_codes[i]) {
      method = i;
      break;
    }
  } 
  
  if (str[0] == 0 || i < 0)  {
    MYWARNING("Ignoring optimization option 'method'. Unrecognized value.\n");
    return;
  }

  for (i = submeth_codes_max[method] - 1; i >= 0; i--)  {
    if (str[1] == submeth_codes[method][i]) {
      submethod = i;
      break;
    }
  } 
  if (str[1] == 0 || i < 0)  {
    MYWARNING("Unrecognized submethod - using default.\n");
  }
}


int OptimizationOptions::gslOptimize( NLSFunction *F, gsl_vector* x_vec, 
        gsl_matrix *v, IterationLogger *itLog ) {
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
  
  switch (this->method) {
  case SLRA_OPT_METHOD_LM:
    gsl_blas_ddot(solverlm->f, solverlm->f, &this->fmin);
    gsl_multifit_gradient(solverlm->J, solverlm->f, g);
    gsl_vector_scale(g, 2);
    {
      gsl_vector *g2 = gsl_vector_alloc(g->size);
      F->computeFuncAndGrad(x_vec, NULL, g2);
      gsl_vector_sub(g2, g);
      if (gsl_vector_max(g2) > 1e-10 || gsl_vector_min(g2) < -1e-10) {
        PRINTF("Gradient error, max = %14.10f,  min = %14.10f  ...",
               gsl_vector_max(g2), gsl_vector_min(g2));
        print_vec(g2);
      }
      gsl_vector_free(g2);
    }
    if (itLog != NULL) {
      itLog->reportIteration(0, solverlm->x, this->fmin, g);
    }
    break;
  case SLRA_OPT_METHOD_QN:
    this->fmin = gsl_multimin_fdfminimizer_minimum(solverqn);
    if (itLog != NULL) {
      itLog->reportIteration(0, solverqn->x, this->fmin, solverqn->gradient);
    }
    break;
  case SLRA_OPT_METHOD_NM:
    this->fmin = gsl_multimin_fminimizer_minimum( solvernm );
    if (itLog != NULL) {
      itLog->reportIteration(this->iter, solvernm->x, this->fmin, NULL);
    }
    break;
  }

  while (status_dx == GSL_CONTINUE && 
	 status_grad == GSL_CONTINUE &&
	 status == GSL_SUCCESS &&
	 this->iter < this->maxiter) {
  	if (this->method == SLRA_OPT_METHOD_LM && this->maxx > 0) {
  	  if (gsl_vector_max(solverlm->x) > this->maxx || 
  	      gsl_vector_min(solverlm->x) < -this->maxx ){
  	    break;
	    }
	  }

    this->iter++;
    switch (this->method) {
    case SLRA_OPT_METHOD_LM: /* Levenberg-Marquardt */
      status = gsl_multifit_fdfsolver_iterate(solverlm);
      gsl_multifit_gradient(solverlm->J, solverlm->f, g);
      gsl_vector_scale(g, 2);

      /* check the convergence criteria */
      if (this->epsabs != 0 || this->epsrel != 0) {
        status_dx = gsl_multifit_test_delta(solverlm->dx, solverlm->x, 
	  				  this->epsabs, this->epsrel);
	  	} else {
	  	  status_dx = GSL_CONTINUE;
	  	}
      status_grad = gsl_multifit_test_gradient(g, this->epsgrad);
      gsl_blas_ddot(solverlm->f, solverlm->f, &this->fmin);
      if (itLog != NULL) {
        itLog->reportIteration(this->iter, solverlm->x, this->fmin, g);
      }
      break;
    case SLRA_OPT_METHOD_QN:
      status = gsl_multimin_fdfminimizer_iterate( solverqn );

      /* check the convergence criteria */
      status_grad = gsl_multimin_test_gradient(
          gsl_multimin_fdfminimizer_gradient(solverqn), this->epsgrad);
      status_dx = gsl_multifit_test_delta(solverqn->dx, solverqn->x, 
	 				 this->epsabs, this->epsrel);  		    
      this->fmin = gsl_multimin_fdfminimizer_minimum(solverqn);      
      if (itLog != NULL) {
        itLog->reportIteration(this->iter, solverqn->x, this->fmin, solverqn->gradient);
      }
      break;
    case SLRA_OPT_METHOD_NM:
      status = gsl_multimin_fminimizer_iterate( solvernm );
      /* check the convergence criteria */
      size = gsl_multimin_fminimizer_size( solvernm );
      status_dx = gsl_multimin_test_size( size, this->epsx );
      this->fmin = gsl_multimin_fminimizer_minimum( solvernm );
      if (itLog != NULL) {
        itLog->reportIteration(this->iter, solvernm->x, this->fmin, NULL);
      }
      break;
    }
  } 
  if (this->iter >= this->maxiter) {
    status = EITER;
  }

  switch (this->method) {
  case  SLRA_OPT_METHOD_LM:
    gsl_vector_memcpy(x_vec, solverlm->x);
    if (v != NULL) {
      gsl_multifit_covar(solverlm->J, this->epscov, v); /* ??? Different eps */
    }
    gsl_blas_ddot(solverlm->f, solverlm->f, &this->fmin);
    break;
  case SLRA_OPT_METHOD_QN:
    gsl_vector_memcpy(x_vec, solverqn->x);
    this->fmin = solverqn->f;
    break;
  case SLRA_OPT_METHOD_NM:
    gsl_vector_memcpy(x_vec, solvernm->x);
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

static void normalizeJacobian( gsl_matrix *jac, gsl_vector *scaling ) {
  for (int i = 0; i < jac->size2; i++) {
    gsl_vector jac_col = gsl_matrix_column(jac, i).vector;
    gsl_vector_set(scaling, i, 1 / gsl_blas_dnrm2(&jac_col));
    gsl_vector_scale(&jac_col, gsl_vector_get(scaling, i));
  }
}			   

static void moveGN( const gsl_matrix *Vt, const gsl_vector *sig2,  const gsl_vector *ufuncsig,
               double lambda2, gsl_vector * dx, int k, gsl_vector * scaling ) {
  gsl_vector_set_zero(dx);
  double threshold = gsl_vector_get(sig2, 0) * Vt->size2 * DBL_EPSILON * 100;
  
  size_t i;
  
  for (i = 0; (i<sig2->size) &&(gsl_vector_get(sig2, i) >= threshold); i++) {
    gsl_vector VtRow = gsl_matrix_const_row(Vt, i).vector;
    gsl_blas_daxpy(gsl_vector_get(ufuncsig, i) /
                   (gsl_vector_get(sig2, i) * gsl_vector_get(sig2, i) + lambda2), &VtRow, dx);
  }
  
/*  if (i >= k) {
    PRINTF("Pseudoinverse threshold exceeded.\n");
    PRINTF("Threshold: %g, i: %d, last: %g, next: %g\n",
           threshold, i, gsl_vector_get(sig2, i-1), gsl_vector_get(sig2, i));
 }*/
 

  
  if (scaling != NULL) {
    gsl_vector_mul(dx, scaling);
  }
}			   

int OptimizationOptions::lmpinvOptimize( NLSFunction *F, gsl_vector* x_vec, 
        IterationLogger *itLog ) {
  int status, status_dx, status_grad, k;
  double g_norm, x_norm;

  if (this->maxiter < 0 || this->maxiter > 5000) {
    throw new Exception("opt.maxiter should be in [0;5000].\n");   
  }
  int scaled = 1; //this->submethod;
  
  /* LM */
  gsl_matrix *jac = gsl_matrix_alloc(F->getNsq(), F->getNvar());
  gsl_vector *func = gsl_vector_alloc(F->getNsq());
  gsl_vector *g = gsl_vector_alloc(F->getNvar());
  gsl_vector *x_cur = gsl_vector_alloc(F->getNvar());
  gsl_vector *x_new = gsl_vector_alloc(F->getNvar());
  gsl_vector *dx = gsl_vector_alloc(F->getNvar());
  gsl_vector *scaling = scaled ? gsl_vector_alloc(F->getNvar()) : NULL;

  gsl_matrix *tempv = gsl_matrix_alloc(jac->size2, jac->size2);
  gsl_vector *tempufuncsig = gsl_vector_alloc(jac->size2);
  gsl_vector *templm = gsl_vector_alloc(jac->size2);
  gsl_vector *sig = gsl_vector_alloc(mymin(jac->size1, jac->size2));

  double lambda2 = 0, f_new;
  int start_lm = 1;
  
  /* Determine optimal work */
  size_t status_svd = 0, minus1 = -1;
  double tmp;
  dgesvd_("A", "O", &jac->size2, &jac->size1, jac->data, &jac->tda, sig->data,
     tempv->data, &tempv->size2, NULL, &jac->size1, &tmp, &minus1, &status_svd);
  gsl_vector *work_vec = gsl_vector_alloc(tmp);
  
  /* optimization loop */
  Log::lprintf(Log::LOG_LEVEL_FINAL, "SLRA optimization:\n");
    
  status = GSL_SUCCESS;  
  status_dx = GSL_CONTINUE;
  status_grad = GSL_CONTINUE;  
  this->iter = 0;
  
  gsl_vector_memcpy(x_cur, x_vec);
  
  F->computeFuncAndJac(x_cur, func, jac);
  gsl_multifit_gradient(jac, func, g);
  gsl_vector_scale(g, 2);
  gsl_blas_ddot(func, func, &this->fmin);
  if (itLog != NULL) {
    itLog->reportIteration(0, x_cur, this->fmin, g);
  }
  
  
  {
    gsl_vector *g2 = gsl_vector_alloc(g->size);
    F->computeFuncAndGrad(x_vec, NULL, g2);
    gsl_vector_sub(g2, g);
    if (gsl_vector_max(g2) > 1e-10 || gsl_vector_min(g2) < -1e-10) {
      PRINTF("Gradient error, max = %14.10f,  min = %14.10f  ...",
             gsl_vector_max(g2), gsl_vector_min(g2));
      print_vec(g2);
    }
    gsl_vector_free(g2);
  }
  
  
  while (status_dx == GSL_CONTINUE &&
         status_grad == GSL_CONTINUE &&
         status == GSL_SUCCESS &&
         this->iter < this->maxiter) {
	/* Check convergence criteria (except dx) */
    if (this->maxx > 0) {
  	  if (gsl_vector_max(x_cur) > this->maxx || gsl_vector_min(x_cur) < -this->maxx ){
  	    break;
      }
    }
  
    this->iter++;

    if (scaling != NULL) {
      normalizeJacobian(jac, scaling);
    }

    
    /* Compute the SVD */
    dgesvd_("A", "O", &jac->size2, &jac->size1, jac->data, &jac->tda, sig->data,
        tempv->data, &tempv->size2, NULL, &jac->size1, work_vec->data, 
		&work_vec->size, &status_svd);

    gsl_blas_dgemv(CblasTrans, -1.0, jac, func, 0.0, tempufuncsig);
    gsl_vector_mul(tempufuncsig, sig);
    while (1) {
      moveGN(tempv, sig, tempufuncsig, lambda2, dx, F->getNEssVar(), scaling);
      gsl_vector_memcpy(x_new, x_cur);
      gsl_vector_add(x_new, dx);
      F->computeFuncAndGrad(x_new, &f_new, NULL);
	  
	    if (f_new <= this->fmin + 1e-16) {
        lambda2 = 0.4 * lambda2;
	      break;
	    }
	    /* Else: update lambda */
	    if (start_lm) {
        lambda2 = gsl_vector_get(sig, 0) * gsl_vector_get(sig, 0);
	      start_lm = 0;
      } else {
        lambda2 = 10 * lambda2;
        Log::lprintf(Log::LOG_LEVEL_ITER, "lambda: %f\n", lambda2);
      }
    }
    /* check the dx convergence criteria */
    if (this->epsabs != 0 || this->epsrel != 0) {
      status_dx = gsl_multifit_test_delta(dx, x_cur, this->epsabs, this->epsrel);
    }     
    gsl_vector_memcpy(x_cur, x_new);

    F->computeFuncAndJac(x_cur, func, jac);
    gsl_multifit_gradient(jac, func, g);
    gsl_vector_scale(g, 2);
    gsl_blas_ddot(func, func, &this->fmin);

    if (itLog != NULL) {
      itLog->reportIteration(this->iter, x_cur, this->fmin, g);
    }
    status_grad = gsl_multifit_test_gradient(g, this->epsgrad);
  } 
  if (this->iter >= this->maxiter) {
    status = EITER;
  }

  gsl_blas_ddot(func, func, &this->fmin);
  
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

  gsl_vector_memcpy(x_vec, x_cur);

  gsl_vector_free(work_vec);
  gsl_matrix_free(jac);
  gsl_vector_free(func);
  gsl_vector_free(g);
  gsl_vector_free(x_cur);
  gsl_vector_free(x_new);
  if (scaling != NULL) {
    gsl_vector_free(scaling);
  }
  gsl_vector_free(dx);
  gsl_matrix_free(tempv);
  gsl_vector_free(tempufuncsig);
  gsl_vector_free(templm);
  gsl_vector_free(sig);
  
  return GSL_SUCCESS; /* <- correct with status */
}


