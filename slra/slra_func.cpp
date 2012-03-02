#include <memory.h>
#include <time.h>

extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}

#include "slra.h"

void allocate_and_prepare_data_reshaped( 
     gsl_matrix* c, int n, 
     opt_and_info *opt, slra_opt_data_reshaped *P, 
     gsl_matrix *perm ) {
  P->m = c->size1;
  P->n = n;
  P->d = c->size2 - n;                            
  P->perm = perm;
                                   
  P->myW = new slraFlexComputationsParams(P->myStruct);
                                   
  /* set other parameters */ \
  P->c = gsl_matrix_alloc(c->size1, c->size2); \
  gsl_matrix_memcpy(P->c, c); \
  P->reggamma = opt->reggamma; \
  /* find Wk  */ \
  P->k = P->myStruct->getK();  \
  P->n_plus_d = P->n + P->d;   \
  P->n_times_d = P->n * P->d;   \
  P->k_times_d = P->myStruct->getK() * P->d; \
  P->k_times_d_times_s = P->k_times_d * P->myW->getS();\
  P->k_times_d_times_s_minus_1 = P->k_times_d_times_s - 1;\
  P->m_times_d = P->m * P->d;\
  P->m_div_k = (int) P->m / P->myStruct->getK();\
  P->s_minus_1 = P->myW->getS() - 1;\
  P->one = 1;\
  P->use_slicot = opt->use_slicot; \
  P->chol_count = 0; \
  P->chol_time = 0; \
//  PREPARE_COMMON_PARAMS(c, n, s, opt, P, 0); 
 
  P->d_times_s = P->d * P->myW->getS();
  P->d_times_m_div_k = P->d * (int) P->m / P->myStruct->getK();
  P->d_times_s_minus_1 = P->d_times_s - 1;

  /* Form reshaped a and b matrices */
//  P->brg_c =  gsl_matrix_alloc(c->size1, c->size2);
  P->brg_c =  c;
  
/*
  gsl_matrix_view src_c = gsl_matrix_view_array(c->data, 
						P->m_div_k, 
						s->k * c->tda);
  gsl_matrix_view src_submat, brg_submat;
  int l;
  for (l = 0; l < s->k; l++) {
    src_submat = gsl_matrix_submatrix(&src_c.matrix, 0, 
				      l * c->tda, P->m_div_k, P->n_plus_d);
    brg_submat = gsl_matrix_submatrix(P->brg_c, P->m_div_k * l, 
				      0, P->m_div_k, P->n_plus_d);
    gsl_matrix_memcpy(&brg_submat.matrix, &src_submat.matrix);
  }*/
//  gsl_matrix_memcpy(P->brg_c, c);


  /* New CholGam */
  P->bx_ext =  gsl_matrix_alloc(P->n_plus_d, P->d);
  P->rb2 = (double*) malloc(P->m_times_d * P->k_times_d_times_s * 
			    sizeof(double));
  P->brg_rb = (double*) malloc(P->m_div_k * P->d * P->d_times_s * 
			       sizeof(double));
  


  P->brg_ldwork = 1 + (P->s_minus_1 + 1)* P->d * P->d +  /* pDW */ 
    3 * P->d + /* 3 * K */
    /* Space needed for MB02CV */
    mymax(P->s_minus_1 + 1, P->m_div_k - 1 - P->s_minus_1) * P->d * P->d;

  /* Modified for comparison */
  P->brg_gamma_vec = (double*) malloc(P->d * (P->d_times_s + P->d) * sizeof(double)); 
  P->brg_gamma = gsl_matrix_alloc(P->d, P->d_times_s + P->d);
  P->brg_dgamma = gsl_matrix_alloc(P->d, P->d_times_s);
  P->brg_tdgamma = gsl_matrix_alloc(P->brg_dgamma->size1, 
				    2 * P->brg_dgamma->size2 - 
				    P->brg_dgamma->size1);
  
  P->brg_dwork  = (double*) malloc((size_t)P->brg_ldwork * sizeof(double));
  P->brg_tmp   = gsl_matrix_alloc(P->d, P->n_plus_d);

  P->brg_yr = gsl_vector_alloc(P->m_times_d);  
  P->brg_f = gsl_vector_alloc(P->m_times_d);  
 
  P->brg_j1b_vec = (double*) calloc(P->d_times_m_div_k * P->d, sizeof(double));
  P->brg_j1b = gsl_matrix_calloc(P->d_times_m_div_k, P->d);

  P->brg_j1b_2 = gsl_matrix_calloc(P->d, P->m * P->d);

  
  P->brg_j1_cvec = gsl_vector_alloc(P->m);
  P->brg_j2_pvec = gsl_vector_alloc(P->n_plus_d);
  
  P->brg_st   = gsl_matrix_alloc(P->m_times_d, P->n_times_d);
  P->brg_jres2  = (double*) malloc( P->m_times_d * sizeof(double));

  P->brg_grad_N_k = gsl_matrix_alloc(P->d, P->d);
  P->brg_grad_Vx_k = gsl_matrix_alloc(P->n, P->d);
  P->brg_grad_tmp1 = gsl_matrix_alloc(P->n, P->d);
  P->brg_grad_tmp2 = gsl_matrix_alloc(P->n, P->d);

  P->brg_perm_tmp = gsl_matrix_alloc(P->n_plus_d, P->d);

 
  
  P->myGamma = new slraFlexGammaComputations(P->myStruct, P->n, opt->use_slicot, opt->reggamma, P->myW);
  P->myDerivative = new slraFlexDerivativeComputations(P->myStruct, P->n, P->myW);
  /*  PRINTF("%p %p %p", P->brg_rb, P->brg_dwork, P->brg_gamma_vec);*/
}

void free_memory_reshaped( slra_opt_data_reshaped *P ) {
  int k;

  gsl_matrix_free(P->c);

    
  gsl_matrix_free(P->brg_c);

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
  gsl_matrix_free(P->brg_j1b_2);

  gsl_vector_free(P->brg_j1_cvec);
  gsl_vector_free(P->brg_j2_pvec);

  free(P->brg_jres2);
  gsl_matrix_free(P->brg_st);
  
  gsl_matrix_free(P->brg_grad_N_k);
  gsl_matrix_free(P->brg_grad_Vx_k);
  gsl_matrix_free(P->brg_grad_tmp1);
  gsl_matrix_free(P->brg_grad_tmp2);

  gsl_matrix_free(P->brg_perm_tmp);
  
  delete P->myStruct;
  delete P->myGamma;
  delete P->myDerivative;
  delete P->myW;
  
  delete P->myCostFun;
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


/* Compute x_ext into params */
static void compute_bxext(const gsl_vector* x, slra_opt_data_reshaped *P) {
  /* Form x_ext */
  xmat2_block_of_xext(gsl_matrix_const_view_vector(x, N, D), P->bx_ext, P->perm, P->brg_perm_tmp);
}

static void compute_reshaped_r(gsl_vector* f, slra_opt_data_reshaped * P) {
  gsl_matrix_view f_mat = gsl_matrix_view_vector(f, M, D); 

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		 P->brg_c, P->bx_ext, 0.0, &f_mat.matrix);
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
/*void reshape_r( gsl_vector *reshaped, gsl_vector * original, void* params, int forward ) {
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

/* ---- Auxiliary functions ---- */

/* 
 *  XMAT2XEXT: x_mat |-> x_ext
 *
 *  x_mat  - the parameter x viewed as a matrix, 
 *  x_ext  = kron( I_k, [ x; -I_d ] ),
 *  params - used for the dimensions n = rowdim(X), 
 *           d = coldim(X), k, and n_plus_d = n + d
 */

void xmat2_block_of_xext(gsl_matrix_const_view x_mat, 
			 gsl_matrix *bx_ext, 
			 gsl_matrix *perm, 
			 gsl_matrix *tmp )
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
  
  if (perm != NULL && tmp != NULL) {
    gsl_matrix_memcpy(tmp, bx_ext);
    
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, 
		   perm, tmp, 0.0, bx_ext);
  }
}


void jacobian_reshaped( slra_opt_data_reshaped* P, gsl_matrix* deriv )
{
  int i, j, l, k, info;
  gsl_matrix_view submat, mat, source; 
  gsl_vector_view subvec, res_vec;
  gsl_matrix_view submat1, submat2, perm_sub_matr;
  gsl_vector_view perm_col;

  /* first term of the Jacobian Gamma^{-1/2} kron(a,I_d) */
/*  for (j = 0; j < N; j++) {
    perm_col = gsl_matrix_column(P->perm, j);
    gsl_blas_dgemv(CblasNoTrans, 1.0, P->brg_c, &perm_col.vector, 0.0, 
        P->brg_j1_cvec);
  
    gsl_matrix_set_zero(P->brg_j1b_2);
    
    for (i = 0; i < M; i++) {
      for (k = 0; k < D; k++) {
        gsl_matrix_set(P->brg_j1b_2, k, k + i*D, 
            gsl_vector_get(P->brg_j1_cvec, i));
      }
    }
    P->myGamma->multiplyInvCholeskyArray(P->brg_j1b_2->data, 1, P->d);

    submat1 = gsl_matrix_submatrix(deriv, 0, j*D, P->m_times_d, D); 
    gsl_matrix_vec_inv(&submat1.matrix, P->brg_j1b_2->data);
  }*/

  gsl_matrix_view st_trans_mat = gsl_matrix_view_array(P->brg_st->data, N * D, M * D);

  gsl_matrix_set_zero(&st_trans_mat.matrix);
  for (j = 0; j < N; j++) {
    perm_col = gsl_matrix_column(P->perm, j);
    gsl_blas_dgemv(CblasNoTrans, 1.0, P->brg_c, &perm_col.vector, 0.0, P->brg_j1_cvec);

    for (i = 0; i < M; i++) {
      for (k = 0; k < D; k++) {
        gsl_matrix_set(&st_trans_mat.matrix, k + j*D, k + i*D, gsl_vector_get(P->brg_j1_cvec, i));
      }
    }
  }
  P->myGamma->multiplyInvCholeskyArray(P->brg_st->data, 1, N * D);
  gsl_matrix_vec_inv(deriv, P->brg_st->data);

/*  for (j = 0; j < N; j++) {

    perm_col = gsl_matrix_column(P->perm, j);
    gsl_blas_dgemv(CblasNoTrans, 1.0, P->brg_c, &perm_col.vector, 0.0, 
        P->brg_j1_cvec);
  
    gsl_matrix_set_zero(P->brg_j1b_2);
    
    for (i = 0; i < M; i++) {
      for (k = 0; k < D; k++) {
        gsl_matrix_set(P->brg_j1b_2, k, k + i*D, 
            gsl_vector_get(P->brg_j1_cvec, i));
      }
    }
    P->myGamma->multiplyInvCholeskyArray(P->brg_j1b_2->data, 1, P->d);

    submat1 = gsl_matrix_submatrix(deriv, 0, j*D, P->m_times_d, D); 
    gsl_matrix_vec_inv(&submat1.matrix, P->brg_j1b_2->data);
  }*/


  /* second term (naive implementation) */
  res_vec = gsl_vector_view_array(P->brg_jres2, P->m_times_d);

  gsl_matrix *bx_ext = P->bx_ext;
  for (i = 0; i < N; i++) {
    for (j = 0; j < D; j++) {
      P->myDerivative->computeDijGammaYr(&res_vec.vector, P->bx_ext, P->perm, 
             i, j, P->brg_yr);
      
      /* solve st_ij = Gamma^{-1/2}st_ij */
      P->myGamma->multiplyInvCholeskyVector(&res_vec.vector, 1);

      subvec = gsl_matrix_column(P->brg_st, i*D+j);
      /* New (nonreshaped) */
      gsl_vector_memcpy(&subvec.vector, &res_vec.vector);
    }
  }

  /* deriv = deriv - 0.5 * st */
  gsl_matrix_scale(P->brg_st, 0.5);
  gsl_matrix_sub(deriv, P->brg_st);
}

void grad_reshaped( slra_opt_data_reshaped* P, gsl_vector* grad )
{
  gsl_matrix_view yr_matr, grad_matr, perm_sub_matr;

  yr_matr = gsl_matrix_view_vector(P->brg_yr, P->m, P->d);
  grad_matr = gsl_matrix_view_vector(grad, P->n, P->d);

  /* Compute term 1 */
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 2.0, P->brg_c, 
                 &yr_matr.matrix, 0.0, P->brg_perm_tmp);
  perm_sub_matr = gsl_matrix_submatrix(P->perm, 0, 0, P->n_plus_d, P->n);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &perm_sub_matr.matrix, 
                 P->brg_perm_tmp, 0.0, &grad_matr.matrix);
  
  /* Compute term 2 */  
  P->myDerivative->computeYrtDgammaYr(P->brg_perm_tmp, P->bx_ext, P->brg_yr);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &perm_sub_matr.matrix, 
                 P->brg_perm_tmp, 0.0, P->brg_grad_tmp1);
  
  gsl_matrix_sub(&grad_matr.matrix, P->brg_grad_tmp1);
}


/* 
 *  slra_F_: STLS cost function evaluation for QN 
 *
 *  x      - row-wise vectorized matrix X
 *  params - parameters for the optimization
 */

double slra_f_reshaped_ (const gsl_vector* x, void* params)
{
  slra_opt_data_reshaped *P = (slra_opt_data_reshaped *) params;
  double ftf;

  /* Use yr as a temporary variable */
  slra_f_reshaped(x, P, P->brg_yr);
  gsl_blas_ddot(P->brg_yr, P->brg_yr, &ftf);

  return ftf;
}

int slra_f_reshaped (const gsl_vector* x, void* params, gsl_vector* f)
{
  slra_opt_data_reshaped *P = (slra_opt_data_reshaped *) params;

  xmat2_block_of_xext(gsl_matrix_const_view_vector(x, N, D), P->bx_ext, P->perm, P->brg_perm_tmp);
  P->myGamma->computeCholeskyOfGamma(P->bx_ext);
  compute_reshaped_r(f, P);
  P->myGamma->multiplyInvCholeskyVector(f, 1);

  return GSL_SUCCESS;
}

int slra_df_reshaped (const gsl_vector* x, void* params, gsl_matrix* deriv)
{
  slra_opt_data_reshaped *P =  (slra_opt_data_reshaped *) params;

  xmat2_block_of_xext(gsl_matrix_const_view_vector(x, N, D), P->bx_ext, P->perm, P->brg_perm_tmp);
  P->myGamma->computeCholeskyOfGamma(P->bx_ext);
  compute_reshaped_r(P->brg_yr, P);
  P->myGamma->multiplyInvGammaVector(P->brg_yr);
  jacobian_reshaped(P, deriv);

  return GSL_SUCCESS;
}


int slra_fdf_reshaped (const gsl_vector* x, void* params, gsl_vector* f, 
		       gsl_matrix* deriv)
{
  slra_opt_data_reshaped *P = (slra_opt_data_reshaped *) params;

  xmat2_block_of_xext(gsl_matrix_const_view_vector(x, N, D), P->bx_ext, P->perm, P->brg_perm_tmp);
  P->myGamma->computeCholeskyOfGamma(P->bx_ext);
  
  compute_reshaped_r(P->brg_yr, P);
  P->myGamma->multiplyInvCholeskyVector(P->brg_yr, 1);

  gsl_vector_memcpy(f, P->brg_yr);

  P->myGamma->multiplyInvCholeskyVector(P->brg_yr, 0);
  
  jacobian_reshaped(P, deriv);

  return GSL_SUCCESS;
}



void slra_fdf_reshaped_ (const gsl_vector* x, void* params, 
			 double *f, gsl_vector* grad)
{
  slra_opt_data_reshaped *P = (slra_opt_data_reshaped *)params;
  gsl_matrix_const_view x_mat = gsl_matrix_const_view_vector( x, N, D );

  xmat2_block_of_xext(gsl_matrix_const_view_vector(x, N, D), P->bx_ext, P->perm, P->brg_perm_tmp);
  P->myGamma->computeCholeskyOfGamma(P->bx_ext);
  
  compute_reshaped_r(P->brg_yr, P);
  if (f !=  NULL) {
    gsl_vector_memcpy(P->brg_f, P->brg_yr);  
  }
  P->myGamma->multiplyInvGammaVector(P->brg_yr);  
  if (f !=  NULL) {
    gsl_blas_ddot(P->brg_f, P->brg_yr, f);
  }

  grad_reshaped(P, grad);

}

void slra_df_reshaped_ (const gsl_vector* x, void* params,  
			gsl_vector* grad)
{
  slra_fdf_reshaped_(x, params, NULL, grad);
}


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
  gsl_vector_view subv;
  gsl_vector_view subp; 	/* subvector of v and p */

  int TM = tt->size1; 		/* = block size */

  gsl_matrix_view submat, source;

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

  int TM = t->size1; 		/* = block size */
  int TN = (t->size2);		/* = s(block size) */


  /* tt - storage for [t_s-1' ... t_1' t_0 t_1 ... t_s-1]. 
     Should be t->size1* (2 * t->size2 - t->size1). */
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

  tmv_prod_new(tt, s, v, m, p);

  free(tt);
}

/* find  w->s */
int get_bandwidth_from_structure( const data_struct* s ) {
  int l, max_nl = 1;

  for (l = 0; l < s->q; l++) {
    if ((!s->a[l].exact) && s->a[l].blocks_in_row > max_nl) {
      max_nl = s->a[l].blocks_in_row;
    }
  }
  
  return max_nl;
}

/* s2w: finds the covariance matrices w from the data structure */
int s2w(const data_struct* s, w_data* w, int n_plus_d, int blocked )
{
  int k, l, i, offset, imax, sum_nl;
  gsl_matrix *zk;
  gsl_matrix_view wi, zkl;
  char err_msg[70];
  int rep;
  int size_wk;
  int ncol;
 
  w->s = get_bandwidth_from_structure(s);
  
  if (blocked) {
    rep = s->k;
    size_wk = s->k * n_plus_d;
  } else {
    rep = 1;
    size_wk = n_plus_d;
  }

  w->a = (gsl_matrix**) malloc(w->s * sizeof(gsl_matrix *));
  zk   = gsl_matrix_alloc(n_plus_d, n_plus_d);
  /* construct w */
  for (k = 0; k < w->s; k++) { 
    gsl_matrix_set_zero(zk);
    sum_nl = 0;
    for (l = 0; l < s->q; l++) { 
      ncol = s->a[l].blocks_in_row * s->a[l].nb;
      zkl = gsl_matrix_submatrix(zk, sum_nl, sum_nl, ncol, ncol); 

      if ((!s->a[l].exact)) {
	offset = s->a[l].nb * k;
	imax   = ncol - offset;
	for (i = 0; i < imax; i++) {
	  if (s->a[l].toeplitz) {
	    gsl_matrix_set(&zkl.matrix, i, i + offset, 1);
          } else {
	    gsl_matrix_set(&zkl.matrix, i + offset, i, 1);
	  }
	}
      } 
      
      sum_nl += ncol;
    }
    w->a[k] = gsl_matrix_calloc(size_wk, size_wk);
    /* w->a[k] = kron(Ik, zk) */
    for (i = 0; i < rep; i++) {
      /* select the i-th diagonal block in a matrix view */
      wi = gsl_matrix_submatrix( w->a[k], i*n_plus_d, 
				 i*n_plus_d, n_plus_d, n_plus_d );
      gsl_matrix_memcpy(&wi.matrix, zk);
    }
  }
  gsl_matrix_free(zk);

  return GSL_SUCCESS;
}

int check_and_adjust_parameters(data_struct *s, flex_struct_add_info *psi) {
  int l;
  psi->total_cols = 0;
  psi->np_scale = 0;
  psi->np_offset = 0;

  for (l = 0; l < s->q; l++) {
    psi->np_offset += s->k * (s->a[l].blocks_in_row - 1) * s->a[l].nb;
    psi->np_scale += s->a[l].nb;
    psi->total_cols += s->a[l].nb * s->a[l].blocks_in_row;
  }  

  return GSL_SUCCESS;
}

#define GET_L(s,l)      (s->a[l].blocks_in_row)
#define GET_T(s,l,m)    (GET_L(s,l) + (m/s->k) - 1) 
#define GET_PLEN(s,l,m) (GET_T(s,l,m) *(s->k) * s->a[l].nb)

int slra_fill_matrix_from_p(gsl_matrix* c,  data_struct *s, gsl_vector* p) {
  int m = c->size1;
  int sum_np = 0, sum_nl = 0;
  int l,j, L, T, p_len;
  gsl_matrix_view p_matr_chunk, c_chunk, p_matr_chunk_sub, c_chunk_sub;
 
  for (l = 0; l < s->q; l++) {
    L = GET_L(s,l);
    T = GET_T(s,l,m);

    p_matr_chunk = gsl_matrix_view_array(&(p->data[sum_np]), 
					 T * s->k, s->a[l].nb);
    c_chunk = gsl_matrix_submatrix(c, 0, sum_nl, m, 
				   s->a[l].nb * s->a[l].blocks_in_row);
    for (j = 0; j < L; j++) {
      p_matr_chunk_sub = gsl_matrix_submatrix(&p_matr_chunk.matrix, 
					      j * s->k, 0, m, s->a[l].nb);
      c_chunk_sub = gsl_matrix_submatrix(&c_chunk.matrix, 0, 
                    (s->a[l].toeplitz ? (L- j -1) * s->a[l].nb : 
		     j * s->a[l].nb), m, s->a[l].nb);
      gsl_matrix_memcpy(&c_chunk_sub.matrix, &p_matr_chunk_sub.matrix);
    }
    sum_np += GET_PLEN(s, l,m);
    sum_nl += s->a[l].nb * s->a[l].blocks_in_row;
  }

  return GSL_SUCCESS;
}


/* Create correction from x */
int slra_correction_reshaped(gsl_vector* p, 
			     data_struct *s, 
			     void* params, 
			     const gsl_vector* x) {
  slra_opt_data_reshaped *P = (slra_opt_data_reshaped *) params;

  /* Compute r(X) */
  slra_f_reshaped(x, params, P->brg_f);
  
  P->myGamma->multiplyInvCholeskyVector(P->brg_f, 0);
  
  
  P->myGamma->correctVector(p, P->myStruct, P->bx_ext, P->brg_f);
}
