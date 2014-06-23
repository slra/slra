#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <cstdarg>

#include "slra.h"

typedef Structure* pStructure;

Structure *createMosaicStructure( gsl_vector * ml, gsl_vector *nk,
               gsl_vector * wk ) {
  enum { ROW_BLW_MOSAIC = 1, BLW_MOSAIC, ELW_MOSAIC } stype;
  
  if (wk == NULL || wk->size == ml->size) {
    stype = ROW_BLW_MOSAIC;
  } else if (wk->size == ml->size * nk->size ) {
    stype = BLW_MOSAIC;
  } else if (wk->size == compute_np(ml,nk)) {
    stype = ELW_MOSAIC;
  } else {
    throw new Exception("Incorrect weight specification\n");   
  }
  
  pStructure *res = new pStructure[nk->size];
  double *pw = (wk == NULL ? NULL : wk->data);
  for (size_t k = 0; k < nk->size; k++) {
    if (stype == ELW_MOSAIC) {
      res[k] = new HLayeredElWStructure(ml->data, ml->size, nk->data[k], pw);
      if (pw != NULL) {
        pw += res[k]->getNp();
      }
    } else {
      res[k] = new HLayeredBlWStructure(ml->data, ml->size, nk->data[k], pw);
      if (stype == BLW_MOSAIC) {
        pw += ml->size;
      }
    } 
  }
  return new StripedStructure(nk->size, res,  stype == ROW_BLW_MOSAIC);
}

size_t compute_np( gsl_vector* ml, gsl_vector *nk ) {
  size_t np = 0;
  size_t i;
  
  for (i = 0; i < ml->size; i++) {
    np += ((size_t)gsl_vector_get(ml, i) - 1) * nk->size; 
  }
  for (i = 0; i < nk->size; i++) {
    np += ((size_t)gsl_vector_get(nk, i)) * ml->size; 
  }
  return np;
}

size_t compute_n( gsl_vector* ml, size_t np ) {
  for (size_t i = 0; i < ml->size; i++) {
    np -= ((size_t)gsl_vector_get(ml, i) - 1);
  }
  if (np <= 0 || (np % ml->size != 0)) {
    throw new Exception("np not compatible with structure specification");
  }

  return np / ml->size;
}

/* gsl_matrix_vectorize: vectorize column-wise a gsl_matrix */
void gsl_matrix_vectorize(double* v, const gsl_matrix* m)
{
  for (size_t j = 0; j < m->size2; j++) {
    for (size_t i = 0; i < m->size1; i++) {
      v[i+j*m->size1] = gsl_matrix_get(m,i,j); 
    }
  }
}

/* gsl_matrix_vec_inv: gsl_matrix from an array */
void gsl_matrix_vec_inv(gsl_matrix* m, const double* v)
{
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      gsl_matrix_set(m, i, j, v[i + j * m->size1]);
    }
  }
}

/* print matrix */
void print_mat(const gsl_matrix* m)
{
  PRINTF("\n");
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++)
      PRINTF("%16.14f ", gsl_matrix_get(m, i, j));
    PRINTF("\n");
  }
  PRINTF("\n");
}


/* print matrix */
void print_mat_tr(const gsl_matrix* m)
{
  PRINTF("\n");
  for (size_t j = 0; j < m->size2; j++) {
    for (size_t i = 0; i < m->size1; i++) {
      PRINTF("%16.14f ", gsl_matrix_get(m, i, j));
    }
    PRINTF("\n");
  }
  PRINTF("\n");
}


/* print_arr: print array */
void print_arr(const double* a, size_t n)
{
  size_t i;

  PRINTF("\n");
  for (i = 0; i < n; i++)
    PRINTF("%f ",*(a+i));
  PRINTF("\n");
}

void print_vec(const gsl_vector* a)
{
  size_t i;

  PRINTF("\n");
  for (i = 0; i < a->size; i++)
    PRINTF("%f ",gsl_vector_get(a, i));
  PRINTF("\n");
}

int read_mat( gsl_matrix *a, const char * filename ) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, "Error opening file %s\n", filename);
    return 0;
  }
  gsl_matrix_fscanf(file, a);
  fclose(file);
  return 1;
}

int read_vec( gsl_vector *a, const char * filename ) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, "Error opening file %s\n", filename);
    return 0;
  }
  gsl_vector_fscanf(file, a);
  fclose(file);
  return 1;
}

int read_vec_uint( gsl_vector_uint *a, const char * filename ) {
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    Log::lprintf(Log::LOG_LEVEL_NOTIFY, "Error opening file %s\n", filename);
    return 0;
  }
  gsl_vector_uint_fscanf(file, a);
  fclose(file);
  return 1;
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

char *strncpy0( char *dest, const char * src, size_t buf_len ) {
  char *res = strncpy(dest, src, buf_len - 1);
  dest[buf_len - 1] = 0;
  return res;
}












