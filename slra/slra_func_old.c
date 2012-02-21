#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>

#include "slra.h"


void allocate_and_prepare_data_old( gsl_matrix* c, int n, const data_struct* s, opt_and_info *opt, slra_opt_data_old *P ) {
  PREPARE_COMMON_PARAMS(c, n, s, opt, P,1); 

  
  /* Preallocate memory for f and df */
  P->x_ext = gsl_matrix_calloc(P->w.a[0]->size1, P->k_times_d);
  P->yr = gsl_vector_alloc(P->m_times_d);  

  /*CholGam */
  P->rb = (double*) calloc(P->m_times_d * P->k_times_d_times_s, sizeof(double));

  P->ldwork = 1 + (P->s_minus_1 + 1)* P->k_times_d * P->k_times_d +  /* pDW */ 
                       3 * P->k_times_d + /* 3 * K */
                       mymax(P->s_minus_1 + 1,  P->m_div_k - 1 - P->s_minus_1) * P->k_times_d * P->k_times_d; /* Space needed for MB02CV */

  P->dwork  = (double*) malloc((size_t)P->ldwork * sizeof(double));
  P->gamma = gsl_matrix_alloc(P->k_times_d, P->k_times_d_times_s);  
  P->gamma_vec = (double*) malloc(P->k_times_d * P->k_times_d_times_s *sizeof(double));


  /* Jacobian*/
  P->jres1 = malloc(P->m_times_d * P->n_times_d * sizeof(double));
  P->jres2  = malloc( P->m_times_d * sizeof(double));

  P->dgamma = gsl_matrix_alloc(P->k_times_d, P->k_times_d_times_s);
  P->st   = gsl_matrix_alloc(P->m_times_d, P->n_times_d);

  P->tmp   = gsl_matrix_alloc(P->k_times_d, P->w.a[0]->size1);

}


void free_memory_old( slra_opt_data_old *P ) {
  int k;

  gsl_matrix_free(P->c);

  /* free the allocated memory for w */
  for (k = 0; k < P->w.s; k++) 
    gsl_matrix_free(P->w.a[k]);
  free(P->w.a);
  
  /* free preallocated memory for computations */
  gsl_matrix_free(P->tmp);
  gsl_matrix_free(P->gamma);
  free(P->dwork);
  free(P->gamma_vec);
  gsl_vector_free(P->yr);
  free(P->rb);

  gsl_matrix_free(P->x_ext);

  free(P->jres1);
  free(P->jres2);
  gsl_matrix_free(P->dgamma);
  gsl_matrix_free(P->st);
}


#define M (P->m)
#define N (P->n)
#define D (P->d)
#define S (P->w.s)
#define SIZE_W (P->w.a[0]->size1)

/* Compute x_ext into params */
static void compute_xext( const gsl_vector* x, slra_opt_data_old *P ) {
  /* reshape x as an nxd matrix x_mat */
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );

  /* Form x_ext */
  xmat2xext( x_mat, P->x_ext, P->k);
}

/* compute f = vec((ax-b)') */
static void compute_f( gsl_vector* f, const gsl_vector* x, slra_opt_data_old *P ) {

  gsl_matrix_view f_mat = gsl_matrix_view_vector(f, M, D); 
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );
  gsl_matrix_view bx_ext = gsl_matrix_submatrix(P->x_ext, 0, 0, N+D, D);

   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
     P->c, &bx_ext.matrix, 0.0, &f_mat.matrix);
 /* gsl_matrix_memcpy(&f_mat.matrix, P->b);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
     P->a, &x_mat.matrix, -1.0, &f_mat.matrix);*/
}

static void compute_c_minus_1_2_f( gsl_vector* f, int trans, slra_opt_data_old *P ) {
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




static void compute_c_minus_1_f( gsl_vector* f, slra_opt_data_old *P ) {
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


/* 
*  slra_F: STLS cost function evaluation 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*  f      - cost function
*/

int slra_f (const gsl_vector* x, void* params, gsl_vector* f)
{
  slra_opt_data_old *P = params;


  compute_xext(x, P);
  compute_f(f, x, P);
  cholgam(P);
  compute_c_minus_1_2_f(f, 1, P);
 
  return GSL_SUCCESS;
}


/* 
*  slra_F_: STLS cost function evaluation for QN 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*/

double slra_f_ (const gsl_vector* x, void* params)
{
  slra_opt_data_old *P = params;

  double ftf;

  /* Use yr as a temporary variable */
  slra_f(x, P, P->yr);
  gsl_blas_ddot(P->yr, P->yr, &ftf);

  return ftf;
}




/* 
*  slra_FD: STLS first derivative evaluation 
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*  deriv  - first derivative
*/

int slra_df (const gsl_vector* x, 
       void* params, gsl_matrix* deriv)
{
  slra_opt_data_old *P = params;

  compute_xext(x, P);
  cholgam(P);
  compute_f(P->yr, x, P);

  compute_c_minus_1_f(P->yr, P);
  jacobian(params, deriv);

  return GSL_SUCCESS;
}


/* 
*  slra_FDF: STLS cost function and first derivative evaluation
*
*  x      - row-wise vectorized matrix X
*  params - parameters for the optimization
*  f      - cost function
*  deriv  - first derivative
*/

int slra_fdf (const gsl_vector* x, void* params, gsl_vector* f, 
        gsl_matrix* deriv)
{
  slra_opt_data_old *P = params;


  compute_xext(x, P);
  compute_f(f, x, P);
  cholgam(P);
  compute_c_minus_1_2_f(f, 1, P);
  gsl_vector_memcpy(P->yr, f);
  compute_c_minus_1_2_f(P->yr, 0, P);
  jacobian(P, deriv);

  return GSL_SUCCESS;
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

void cholgam( slra_opt_data_old* P )
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
       P->x_ext, P->w.a[k], 0.0, P->tmp);
    /* compute tmp * x_ext */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
       P->tmp, P->x_ext, 0.0, &submat.matrix);
  }

  gsl_matrix_vectorize(P->gamma_vec, P->gamma);
  
  
  /*PRINTF("%p %p %p\n", P->rb, P->dwork, P-> gamma_vec);*/
  /*PRINTF("cholgam: PDW = %d, LDWORK = %d, 3K = %d\n, s-1 = %d", 1 + (P->s_minus_1 + 1)* P->k_times_d * P->k_times_d, ldwork, 3 * P->k_times_d, P->s_minus_1 );*/

  /* Cholesky factorization of Gamma */
  mb02gd_("R", "N",  
    &P->k_times_d,   /* block size */
    &P->m_div_k,    /* block_dim(GAMMA) */
    &P->s_minus_1,   /* non-zero block superdiagonals */
    &zero,     /* previous computations */
    &P->m_div_k,           /* to-be-computed */
    P->gamma_vec, /* non-zero part of first block row of GAMMA */
    &P->k_times_d,        /* row_dim(gamma) */
    P->rb,       /* packed Cholesky factor */
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

void jacobian( slra_opt_data_old* P, gsl_matrix* deriv )
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
        gsl_matrix_set(&submat.matrix, l, l, gsl_matrix_get(P->c, i, j));
    }
  
  /* res1 = vec(deriv), FORTRAN column-major convention */
  gsl_matrix_vectorize(P->jres1, deriv);

  /* res =  Gamma^{-1/2} * res  */
  dtbtrs_("U", "T", "N", 
    &P->m_times_d, 
    &P->k_times_d_times_s_minus_1, 
    &P->n_times_d, 
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
          w_k_row = gsl_matrix_row (P->w.a[k], l*P->n_plus_d+i);
          gsl_blas_dgemv(CblasTrans, 1.0, P->x_ext, &w_k_row.vector, 0.0, &tmp1_row.vector); 
        }
        
        /* compute submat = submat  + x_ext' * tmp * dx_ext * /
        /* Iterate over rows of dx_ext' * w_k' */
        for (l = 0; l < P->k; l++) {
          tmp1_col = gsl_matrix_column (&submat.matrix, l*D+j);
          w_k_col = gsl_matrix_column (P->w.a[k], l*P->n_plus_d+i);
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



void xmat2xext( gsl_matrix_const_view x_mat, gsl_matrix *x_ext, int K )
{
  int i;
  gsl_matrix_view submat, source;
  int n = x_mat.matrix.size1,
      d = x_mat.matrix.size2;


  /* set block (1,1) of x_ext to [ x_mat; -I_d ] */
  xmat2_block_of_xext(x_mat, x_ext, NULL, NULL);

  /* copy block (1,1) in (2,2),...,(k,k) */
  source = gsl_matrix_submatrix(x_ext, 0, 0, n+d, d);
  for (i = 1; i < K; i++) {
    submat = gsl_matrix_submatrix(x_ext, i*(n+d), i*d, (n+d), d);
    gsl_matrix_memcpy(&submat.matrix, &source.matrix);
  }
}


