#include "stls.h"

/* ************************************************ */
/* m_to_gsl_matrix: convert the Matlab style column */
/* major mxn matrix array to a GSL matrix           */
/* ************************************************ */

void m_to_gsl_matrix(gsl_matrix* a_gsl, double* a_m) 
{
  int i, j;

  for (i = 0; i < a_gsl->size1; i++) {
    for (j = 0; j < a_gsl->size2; j++) {
      gsl_matrix_set(a_gsl, i, j, a_m[i + j * a_gsl->size1]);
    }
  }
}

/* ************************************************ */
/* gsl_to_m_matrix: convert the GSL mxn matrix to a */
/* Matlab style column major matrix array           */
/* ************************************************ */

void gsl_to_m_matrix(double* a_m, gsl_matrix* a_gsl) 
{
  int i, j;

  for (i = 0; i < a_gsl->size1; i++) {
    for (j = 0; j < a_gsl->size2; j++) {
      a_m[i + j * a_gsl->size1] = gsl_matrix_get(a_gsl, i, j);
    }
  }
}
