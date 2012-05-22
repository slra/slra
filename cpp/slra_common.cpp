#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <cstdarg>

#include "slra.h"


Exception::Exception( const char *format, ... ) { 
  va_list vl;
  va_start(vl, format);  
  myMsg[MSG_MAX-1] = 0;
  vsnprintf(myMsg, MSG_MAX-1, format, vl); 
}

Structure *createMosaicStructure( gsl_vector * ml,  gsl_vector *nk, 
               gsl_vector * wk, int np_comp ) {
  if (wk == NULL || wk->size == ml->size * nk->size || wk->size == ml->size) { 
    return new MosaicHStructure(ml, nk, wk);
  } 
  if (wk->size == np_comp) {
    return new WMosaicHStructure(ml, nk, wk);
  } 
  throw new Exception("Incorrect weight specification\n");   
}

void OptimizationOptions::str2Method( const char *str )  {
  char meth_codes[] = "lqn", 
       sm_codes_lm[] = "ls", sm_codes_qn[] = "b2pf", sm_codes_nm[] = "n2r";
  char *submeth_codes[] = { sm_codes_lm, sm_codes_qn, sm_codes_nm };

  int submeth_codes_max[] = { 
    sizeof(sm_codes_lm) / sizeof(sm_codes_lm[0]) - 1, 
    sizeof(sm_codes_qn) / sizeof(sm_codes_qn[0]) - 1, 
    sizeof(sm_codes_nm) / sizeof(sm_codes_nm[0]) - 1
  };
  int meth_code_max = sizeof(submeth_codes_max) / sizeof(submeth_codes_max[0]);
  int i;
  
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

void OptimizationOptions::str2Disp( const char *str )  {
  char *str_disp[] = {"notify", "final", "iter", "off" };
  int i;
  
  for (i = 0; i < sizeof(str_disp) / sizeof(str_disp[0]); i++) {
    if (strcmp(str_disp[i], str) == 0) {
      disp = i;
    }
  }
    
  if (i < 0) {
    MYWARNING("Ignoring optimization option 'disp'. Unrecognized value.\n");
  }
}


/* gsl_matrix_vectorize: vectorize column-wise a gsl_matrix */
void gsl_matrix_vectorize(double* v, const gsl_matrix* m)
{
  for (int j = 0; j < m->size2; j++) {
    for (int i = 0; i < m->size1; i++) {
      v[i+j*m->size1] = gsl_matrix_get(m,i,j); 
    }
  }
}

/* gsl_matrix_vec_inv: gsl_matrix from an array */


void gsl_matrix_vec_inv(gsl_matrix* m, const double* v)
{
  for (int i = 0; i < m->size1; i++) {
    for (int j = 0; j < m->size2; j++) {
      gsl_matrix_set(m, i, j, v[i + j * m->size1]);
    }
  }
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
void print_arr(const double* a, int n)
{
  int i;

  PRINTF("\n");
  for (i = 0; i < n; i++)
    PRINTF("%f ",*(a+i));
  PRINTF("\n");
}

void print_vec(const gsl_vector* a)
{
  int i;

  PRINTF("\n");
  for (i = 0; i < a->size; i++)
    PRINTF("%f ",gsl_vector_get(a, i));
  PRINTF("\n");
}

int compute_np( gsl_vector* ml, gsl_vector *nk ) {
  int np = 0;
  int i;
  
  for (i = 0; i < ml->size; i++) {
    np += ((int)gsl_vector_get(ml, i) - 1) * nk->size; 
  }
  for (i = 0; i < nk->size; i++) {
    np += ((int)gsl_vector_get(nk, i)) * ml->size; 
  }
  return np;
}


int compute_n( gsl_vector* ml, int np ) {
  for (int i = 0; i < ml->size; i++) {
    np -= ((int)gsl_vector_get(ml, i) - 1);
  }

  if (np <= 0 || (np % ml->size != 0)) {
    throw new Exception("np not compatible with structure specification");
  }

  return np / ml->size;
}

void Cholesky::multInvCholeskyTransMatrix( gsl_matrix * yr_matr, int trans ) { 
  for (int i = 0; i < yr_matr->size1; i++) {
    gsl_vector_view row = gsl_matrix_row(yr_matr, i);
    multInvCholeskyVector(&row.vector, trans);
  }
}

const gsl_vector *vecChkNIL( const gsl_vector &vec ) {
  return  vec.data != NULL ? &vec : NULL; 
}

gsl_vector *vecChkNIL( gsl_vector &vec  ) {
  return  vec.data != NULL ? &vec : NULL; 
}

gsl_matrix *matChkNIL( gsl_matrix &mat_vw ) {
  return  mat_vw.data != NULL ? &mat_vw : NULL; 
}

void tolowerstr( char * str ) {
  char *c;
  for (c = str; *c != '\0'; c++) {
    *c = tolower(*c);
  }
} 

void id_kron_a( const gsl_matrix *A, int d,  gsl_matrix *IkronA ) {
  gsl_matrix_set_zero(IkronA);
  for (int  k = 0; k < d; k++) {
    gsl_matrix subA = gsl_matrix_submatrix(IkronA, k * A->size1, k * A->size2, 
                           A->size1, A->size2).matrix;
    gsl_matrix_memcpy(&subA, A);
  }
}


void ls_solve( const gsl_matrix *A, const gsl_matrix *B, gsl_matrix *X ) {
  size_t d = B->size2, i, j, one = 1, lwork = -1, info;
  gsl_matrix *IkronA = gsl_matrix_alloc(A->size1 * d, A->size2 * d);
  double tmp;
  
  double *vecIkronA = new double[IkronA->size1 * IkronA->size2],
         *vecB = new double[B->size1 * B->size2]; 
                          
  id_kron_a(A, d, IkronA);                                      
  gsl_matrix_vectorize(vecIkronA, IkronA);
  gsl_matrix_vectorize(vecB, B);
  
  dgels_("N", &IkronA->size1, &IkronA->size2, &one, vecIkronA,
         &IkronA->size1, vecB, &B->size1, &tmp, &lwork, &info);
  double *work = new double[lwork = tmp];
  dgels_("N", &IkronA->size1, &IkronA->size2, &one, vecIkronA,
         &IkronA->size1, vecB, &B->size1, work, &lwork, &info);

  delete [] work;
  
  gsl_matrix_free(IkronA);
  gsl_matrix_vec_inv(X, vecB);
  
  delete [] vecB;
  delete [] vecIkronA;
}



