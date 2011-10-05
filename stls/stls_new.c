#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "stls.h"

void allocate_and_prepare_data_new( gsl_matrix* a, gsl_matrix* b, const data_struct* s, stls_opt_data *P ) {
  PREPARE_COMMON_PARAMS(a, b, s, P, 1); 
  
  P->m = a->size1;
  P->n = a->size2;
  P->d = b->size2;
 
  P->d_times_s = P->d * P->w.s;
  P->d_times_m_div_k = P->d* (int) P->m / s->k;
  P->d_times_s_minus_1 = P->d_times_s - 1;


  /* Form reshaped a and b matrices */
  P->brg_a =  gsl_matrix_alloc(a->size1, a->size2);
  P->brg_b =  gsl_matrix_alloc(b->size1, b->size2);
  
  gsl_matrix_view src_a = gsl_matrix_view_array(a->data, P->m_div_k, s->k * P->n);
  gsl_matrix_view src_b = gsl_matrix_view_array(b->data, P->m_div_k, s->k * P->d);
  gsl_matrix_view src_submat, brg_submat;
  int l;
  for (l = 0; l < s->k; l++) {
  	src_submat = gsl_matrix_submatrix(&src_a.matrix, 0, l * P->n, P->m_div_k, P->n);
  	brg_submat = gsl_matrix_submatrix(P->brg_a, P->m_div_k * l, 0, P->m_div_k, P->n);
  	gsl_matrix_memcpy(&brg_submat.matrix, &src_submat.matrix);

  	src_submat = gsl_matrix_submatrix(&src_b.matrix, 0, l * P->d, P->m_div_k, P->d);
  	brg_submat = gsl_matrix_submatrix(P->brg_b, P->m_div_k * l, 0, P->m_div_k, P->d);
  	gsl_matrix_memcpy(&brg_submat.matrix, &src_submat.matrix);
  }
 
  
  

  /* New CholGam */
  P->bx_ext =  gsl_matrix_alloc(P->n_plus_d, P->d);

  P->rb2 = (double*) malloc(P->m_times_d * P->k_times_d_times_s * sizeof(double));
  P->brg_rb = (double*) malloc(P->m_div_k * P->d * P->d_times_s * sizeof(double));
  


  P->brg_ldwork = 1 + (P->s_minus_1 + 1)* P->d * P->d +  /* pDW */ 
                       3 * P->d + /* 3 * K */
                       mymax(P->s_minus_1 + 1,  P->m_div_k - 1 - P->s_minus_1) * P->d * P->d; /* Space needed for MB02CV */

  P->brg_gamma_vec = (double*) malloc(P->d * P->d_times_s * sizeof(double));
  P->brg_gamma = gsl_matrix_alloc(P->d, P->d_times_s);
  P->brg_dgamma = gsl_matrix_alloc(P->d, P->d_times_s);
  P->brg_tdgamma = gsl_matrix_alloc(P->brg_dgamma->size1, 2*P->brg_dgamma->size2 - P->brg_dgamma->size1);
  
  P->brg_dwork  = (double*) malloc((size_t)P->brg_ldwork * sizeof(double));
  P->brg_tmp   = gsl_matrix_alloc(P->d, P->n_plus_d);

  P->brg_yr = gsl_vector_alloc(P->m_times_d);  
  P->brg_f = gsl_vector_alloc(P->m_times_d);  

 
  P->brg_j1b_vec = calloc(P->d_times_m_div_k * P->d, sizeof(double));
  P->brg_j1b = gsl_matrix_calloc(P->d_times_m_div_k, P->d);
  
  P->brg_st   = gsl_matrix_alloc(P->m_times_d, P->n_times_d);
  P->brg_jres2  = malloc( P->m_times_d * sizeof(double));



/*  PRINTF("%p %p %p", P->brg_rb, P->brg_dwork, P->brg_gamma_vec);*/
}

void free_memory_new( stls_opt_data *P ) {
  int k;

  for (k = 0; k < P->w.s; k++) 
    gsl_matrix_free(P->w.a[k]);
  free(P->w.a);
  
  
  gsl_matrix_free(P->brg_a);
  gsl_matrix_free(P->brg_b);

  gsl_matrix_free(P->bx_ext);
  gsl_matrix_free(P->brg_tmp);
  gsl_matrix_free(P->brg_gamma);
  gsl_matrix_free(P->brg_dgamma);
  gsl_matrix_free(P->brg_tdgamma);
  
  free(P->brg_dwork);
  free(P->brg_gamma_vec);
  gsl_vector_free(P->brg_yr);
  gsl_vector_free(P->brg_f);
  free(P->rb2);
  free(P->brg_j1b_vec);
  gsl_matrix_free(P->brg_j1b);


  free(P->brg_jres2);
  gsl_matrix_free(P->brg_st);
}

#define M (P->m)
#define N (P->n)
#define D (P->d)

/*
#define M (P->a->size1)
#define N (P->a->size2)
#define D (P->b->size2)*/
#define S (P->w.s)
#define SIZE_W (P->w.a[0]->size1)




/* Compute bx_ext into params */
static void compute_bxext( const gsl_vector* x, stls_opt_data * P ) {
  /* reshape x as an nxd matrix x_mat */
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );

  /* Form x_ext */
  xmat2bxext( x_mat, P->bx_ext );
}




static void compute_reshaped_f( gsl_vector* f, const gsl_vector* x, stls_opt_data * P ) {
  gsl_matrix_view f_mat = gsl_matrix_view_vector(f, M, D); 
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );
 
  gsl_matrix_memcpy(&f_mat.matrix, P->brg_b);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		 P->brg_a, &x_mat.matrix, -1.0, &f_mat.matrix);
}





static void compute_reshaped_c_minus_1_2_f( gsl_vector* f, int trans, stls_opt_data * P ) {
  int info;
  int i;

  dtbtrs_("U", (trans ? "T" : "N"), "N", 
	  &P->d_times_m_div_k, 
	  &P->d_times_s_minus_1, 
	  &P->k, 
	  P->brg_rb, 
	  &P->d_times_s, 
	  f->data, 
	  &P->d_times_m_div_k, 
	  &info);
}

static void compute_reshaped_c_minus_1_f( gsl_vector* f, stls_opt_data * P ) {
  int info;
  
    dpbtrs_("U", 
	    &P->d_times_m_div_k, 
	    &P->d_times_s_minus_1, 
	    &P->k, 
	    P->brg_rb, 
	    &P->d_times_s, 
	    f->data, 
	    &P->d_times_m_div_k, 
	    &info);
  
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
/*void reshape_f( gsl_vector *reshaped, gsl_vector * original, void* params, int forward ) {
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
}*/










int stls_f_new (const gsl_vector* x, void* params, gsl_vector* f)
{
  stls_opt_data *P = params;
  
  compute_bxext(x, P);
  cholbrg(P);
  
  compute_reshaped_f(P->brg_f, x, P);
  compute_reshaped_c_minus_1_2_f(P->brg_f, 1, P);
  gsl_vector_memcpy(f, P->brg_f);

  return GSL_SUCCESS;
}



/* 
*  STLS_F_: STLS cost function evaluation for QN 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*/

double stls_f_new_ (const gsl_vector* x, void* params)
{
  stls_opt_data *P = params;

  double ftf;

  /* Use yr as a temporary variable */
  stls_f_new(x, P, P->brg_yr);
  gsl_blas_ddot(P->brg_yr, P->brg_yr, &ftf);

  return ftf;
}





int stls_df_new (const gsl_vector* x, void* params, gsl_matrix* deriv)
{
  stls_opt_data *P = params;

  compute_bxext(x, params);
  cholbrg(P);
  compute_reshaped_f(P->brg_yr, x, P);
  compute_reshaped_c_minus_1_f(P->brg_yr, P);

  jacobian_new(P, deriv);

  return GSL_SUCCESS;
}






int stls_fdf_new (const gsl_vector* x, void* params, gsl_vector* f, 
	      gsl_matrix* deriv)
{
  stls_opt_data *P = params;

  compute_bxext(x, P);
  cholbrg(P);

  compute_reshaped_f(P->brg_yr, x, P);
  compute_reshaped_c_minus_1_2_f(P->brg_yr, 1, P);

  gsl_vector_memcpy(f, P->brg_yr);

  compute_reshaped_c_minus_1_2_f(P->brg_yr, 0, P);
  
  jacobian_new(P, deriv);

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



void xmat2bxext( gsl_matrix_const_view x_mat, gsl_matrix *bx_ext )
{
  gsl_vector_view diag;
  gsl_matrix_view submat;
  int n = x_mat.matrix.size1,
      d = x_mat.matrix.size2;



  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  submat = gsl_matrix_submatrix(bx_ext, 0, 0, n, d); /* select x in (1,1) */
  gsl_matrix_memcpy(&submat.matrix, &x_mat.matrix); /* assign x */
  submat = gsl_matrix_submatrix(bx_ext, n, 0, d, d); /* select -I in (1,1)*/
  gsl_matrix_set_all(&submat.matrix,0);
  diag   = gsl_matrix_diagonal(&submat.matrix);     /* assign -I */
  
  
  gsl_vector_set_all(&diag.vector, -1);
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

void cholbrg( stls_opt_data* P )
{
  int k, info;
  gsl_matrix_view submat;
  gsl_matrix_view  b_w_k; 
  gsl_matrix *bx_ext = P->bx_ext;


/*  printf("cholbrg: %p\n", P);*/


  const int zero = 0;
  /* MB02GD has very bad description of parameters */


  /* Take block of x_ext - b_x_ext = x_ext(1:(n+d),1:d) */
  
/*  gsl_matrix_view b_x_ext;
  b_x_ext =  gsl_matrix_submatrix(P->x_ext, 0, 0, P->n_plus_d, D);
  bx_ext = &b_x_ext.matrix; */
  
  
/*  PRINTF("x_ext_b:");
  print_mat(bx_ext);*/
  
  

  /* compute brgamma_k = b_x_ext' * w_k(1:(n+d),1:(n+d)) * b_x_ext */
  for (k = 0; k < S; k++) {
    submat = gsl_matrix_submatrix(P->brg_gamma, 0, k*D, D, D);
    b_w_k = gsl_matrix_submatrix(P->w.a[k], 0, 0, P->n_plus_d, P->n_plus_d);
				  
				  
    /* compute tmp = x_ext' * w_k */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, bx_ext, &b_w_k.matrix, 0.0, P->brg_tmp);
    /* compute brgamma_k = tmp * x_ext */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,  P->brg_tmp, bx_ext, 0.0, &submat.matrix);
  }

  gsl_matrix_vectorize(P->brg_gamma_vec, P->brg_gamma);

/*  PRINTF("%p %p %p", P->brg_rb, P->brg_dwork, P->brg_gamma_vec);*/

  
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









void jacobian_new( stls_opt_data* P, gsl_matrix* deriv )
{
  int i, j, l, k, info;
  const int zero = 0, one = 1;
  gsl_matrix_view submat, mat, source; 
  gsl_vector_view subvec, res_vec;

    
  
  
  

	gsl_vector_view tmp1_row, tmp1_col;
	gsl_vector_view w_k_row, w_k_col;


  /* first term of the Jacobian Gamma^{-1/2} kron(a,I_d) */



  gsl_matrix_view submat1, submat2;
  int ibr, ik, i1, j1, kbr;
  double a_kj;
  
  

 /* gsl_matrix_set_zero(P->brg_j1b);*/
  
  
  for (j = 0; j < N; j++) {
    for (ik = 0; ik < P->k; ik++) {

      /* Fill right-hand matrix */
  
      for (ibr = 0; ibr < P->m_div_k; ibr++) {
        for (k = 0; k < D; k++) {
          gsl_matrix_set(P->brg_j1b, k + ibr*D, k, gsl_matrix_get(P->a, ik + ibr * P->k,j)) ;
        }
      }
      

      gsl_matrix_vectorize(P->brg_j1b_vec, P->brg_j1b);
      
      
      dtbtrs_("U", "T", "N", 
	      &P->d_times_m_div_k, 
	      &P->d_times_s_minus_1, 
	      &P->d, 
  	    P->brg_rb, 
	      &P->d_times_s, 
	      P->brg_j1b_vec, 
	      &P->d_times_m_div_k, 
  	    &info);

      submat1 = gsl_matrix_submatrix(deriv, 0 + ik * P->d_times_m_div_k, j*D, P->d_times_m_div_k, D); 

      gsl_matrix_vec_inv(&submat1.matrix, P->brg_j1b_vec);
    }
  }
    

/*  gsl_matrix_memcpy(deriv, P->st); */
/*  PRINTF("jacobian_new: step1 max, min = %f, %f, \n", gsl_matrix_max(P->st), gsl_matrix_min(P->st));*/

  /* second term (naive implementation) */
  res_vec = gsl_vector_view_array(P->brg_jres2, P->m_times_d);

  gsl_matrix *bx_ext = P->bx_ext;
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {

      /* Alternative */

		  /* form dgamma = d/d x_ij gamma */
			gsl_matrix_set_zero(P->brg_tdgamma);


			for (k = 0; k < S; k++) {
			  gsl_matrix_view  b_w_k = gsl_matrix_submatrix(P->w.a[k], 0, 0, P->n_plus_d, P->n_plus_d);
			
	  	  submat = gsl_matrix_submatrix(P->brg_tdgamma, 0, (S-1)*D + k*D, D, D);

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

			for (l = 0; l < S-1; l++) {
				submat = gsl_matrix_submatrix(P->brg_tdgamma, 0, l*D, D, D);
				source = gsl_matrix_submatrix(P->brg_tdgamma, 0, (2* S-2-l)*D, D, D);
				gsl_matrix_transpose_memcpy(&submat.matrix, &source.matrix);
		  }
     	
     	
    	/* compute st_ij = DGamma * yr */
      subvec = gsl_matrix_column(P->brg_st, i*D+j);
    		
      for (l = 0; l < P->k; l++) { 
        gsl_vector_view subvec_res = gsl_vector_subvector(&res_vec.vector, l * P->d_times_m_div_k, P->d_times_m_div_k);
        gsl_vector_view subvec_yr = gsl_vector_subvector(P->brg_yr, l * P->d_times_m_div_k, P->d_times_m_div_k);
				
				
        tmv_prod_new(P->brg_tdgamma, S, &subvec_yr.vector, P->m_div_k, &subvec_res.vector);  
        
    	}
      	/* solve st_ij = Gamma^{-1/2}st_ij */
      	dtbtrs_("U", "T", "N", 
			  		&P->d_times_m_div_k, 
				  	&P->d_times_s_minus_1, 
  					&P->k, 
	  				P->brg_rb, 
		  			&P->d_times_s, 
			  		res_vec.vector.data,
				  	&P->d_times_m_div_k, 
    				&info); 

			/* New (nonreshaped) */
 		  gsl_vector_memcpy(&subvec.vector, &res_vec.vector);
    }
  }




  /* deriv = deriv - 0.5 * st */
  gsl_matrix_scale(P->brg_st, 0.5);
  gsl_matrix_sub(deriv, P->brg_st);

}


