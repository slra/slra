/* test.c: test program for STLS 
   .\test i - test example i, 1 <= i <= 24
   .\test   - test all examples            */ 

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>
#include "stls.h"

/* default constants for the exit condition */
#define MAXITER 600
#define EPSABS  0
#define EPSREL  1e-5
#define EPSGRAD 1e-5
#define DISP    3     /* per iteration */

#define TEST_NUM 7

#define MAX_FN_LEN  60

int read_mat( gsl_matrix *a,  char * filename, FILE * log ) {
  FILE * file;
  file = fopen(filename,"r");
  if (file == NULL ) {
    fprintf(log, "Error opening file %s\n", filename);

    return 0;
  }
  gsl_matrix_fscanf(file,a);
  fclose(file);

  return 1;
}


int read_vec( gsl_vector *a,  char * filename, FILE * log ) {
  FILE * file;
  file = fopen(filename,"r");
  if (file == NULL ) {
    fprintf(log, "Error opening file %s\n", filename);

    return 0;
  }
  gsl_vector_fscanf(file,a);
  fclose(file);

  return 1;
}



void run_test( FILE * log, char * testname, double & time, double & fmin, double &fmin2, int & iter, double &diff, bool silent = false ) {
  gsl_matrix *xt = NULL, *x = NULL, *a = NULL, *b = NULL, *v = NULL;
  gsl_vector *p = NULL, * p2 = NULL;
  
  char str_codes[] = " THUE";
  data_struct s = {9599, 1,1,{'H',16,8}}; /* {1,2,{'T',10,1,'U',1,1}}; */
  opt_and_info opt = {MAXITER, DISP, EPSREL, EPSABS, EPSGRAD, 0.01,  0, 0.0, 0.0};
  int i, j, m = 9599, n = 12, d = 4, tmp, np = 9599;
  int x_given;
  char faname[MAX_FN_LEN], fbname[MAX_FN_LEN], fxname[MAX_FN_LEN], fpname[MAX_FN_LEN],
       fxtname[MAX_FN_LEN], fsname[MAX_FN_LEN], fxresname[MAX_FN_LEN],
       fpresname[MAX_FN_LEN];
  FILE *file;

  /* default file names */
  sprintf(faname,"a%s.txt",testname);
  sprintf(fbname,"b%s.txt",testname);
  sprintf(fxname,"x%s.txt",testname);
  sprintf(fpname,"p%s.txt",testname);
  sprintf(fxtname,"xt%s.txt",testname);
  sprintf(fsname,"s%s.txt",testname);
  sprintf(fxresname,"res_x%s.txt",testname);
  sprintf(fpresname,"res_p%s.txt",testname);

  try {
    if (!silent) {
      fprintf(log, "Running test %s\n", testname);  
    }
   
    /* Read structure */
    file = fopen(fsname,"r");
    if (file == NULL) {
      fprintf(log, "Error opening file %s\n", fsname);
      throw 1;
    }
    fscanf(file, "%d %d %d %d %d %d", &s.m, &n, &d, &s.k, &s.q, &np); 
    m = s.m;
    if ((s.k <= 0) || (m % s.k != 0)) {
      fprintf(log, "Bad k: %d \n", s.k);
    }

    if ((s.q <= 0) || (s.q > 10)) {
	    fprintf(log, "Bad k: %d \n", s.k);
	    throw 1;
    }
    
    for (i = 0; i < s.q; i++)  {
      fscanf(file, "%d %d %d", &tmp, &(s.a[i].ncol), &(s.a[i].nb)); 
      s.a[i].type = str_codes[tmp]; 
    }
    fclose(file);

    /* read the data a, b from a.txt and b.txt */
    if (((a = gsl_matrix_alloc(m, n)) == NULL) || !read_mat(a, faname, log)) {
      throw 1;
    }
    if (((b = gsl_matrix_alloc(m, d)) == NULL) || !read_mat(b, fbname, log)) {
      throw 1;
    }


    if (((x = gsl_matrix_calloc(n, d)) == NULL)) {
      throw 1;
    }
    
    x_given = read_mat(x, fxname, log);


    if (((p = gsl_vector_alloc(np)) == NULL)) { 
      throw 1;
    }

    if (((p2 = gsl_vector_alloc(np)) == NULL)) { 
      throw 1;
    }

    
    if (!read_vec(p, fpname, log)) {
      gsl_vector_free(p);
      p = NULL;
    }
    
    gsl_vector_memcpy(p2,p);

    
//    print_mat(x);

    if (((xt = gsl_matrix_calloc(n, d)) == NULL)) {
      throw 1;
    }
    read_mat(xt, fxtname, log);

    if (((v = gsl_matrix_alloc(n*d, n*d)) == NULL)) {
      throw 1;
    }

    if (silent) {
      opt.disp = 0;
    }


    

    /* call stls */  
    stls(a, b, &s, x, v, &opt, p, x_given, 1 /* Compute correction */);

    if (!silent) {
      print_mat(x);
    }


    file = fopen(fpresname,"w");
    gsl_vector_fprintf(file, p, "%.10f");
    fclose(file);


    gsl_vector_sub(p, p2);
    double dp_norm = gsl_blas_dnrm2(p);
    
  

    file = fopen(fxresname,"w");
    gsl_matrix_fprintf(file, x, "%.10f");
    fclose(file);
    

    diff = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j < d; j++) {
        double td;  
        td = gsl_matrix_get(xt, i, j) - gsl_matrix_get(x, i, j);
        diff += td*td;
      }
    }
    diff=sqrt(diff);
    
    /* print and save result */
    if (!silent) {
      fprintf(log, "Result: (time, fmin, diff) = %10.8f %10.8f %f\n", opt.time, opt.fmin, diff);
    }
    time = opt.time;
    fmin = opt.fmin;
    iter = opt.iter;
    fmin2 = dp_norm * dp_norm;
    

    
    throw 0;
  } catch(...) {
 
    /* free the memory */
    if (a != NULL) {
      gsl_matrix_free(a);
    }
    if (b != NULL) {
      gsl_matrix_free(b);
    }
    if (x != NULL) {
      gsl_matrix_free(x);
    }
    if (xt != NULL) {
      gsl_matrix_free(xt);
    }
    if (v != NULL) {
      gsl_matrix_free(v);
    }
    if (p != NULL) {
      gsl_vector_free(p);
    }
    if (p2 != NULL) {
      gsl_vector_free(p2);
    }
  }
}




int main(int argc, char *argv[])
{
  double times[TEST_NUM+1], misfits[TEST_NUM+1], misfits2[TEST_NUM+1],  diffs[TEST_NUM+1];
  int iters[TEST_NUM+1];
  char num[10];
  int  i;

  if (argc == 2) {
    run_test(stdout, argv[1], times[0], misfits[0], misfits2[0], iters[0], diffs[0]);
  } else { /* test all examples */
    printf("\n------------------ Testing all examples  ------------------\n\n");

    for( i = 1; i <= TEST_NUM; i++ ) {
      sprintf(num, "%d", i);
      run_test(stdout, num, times[i], misfits[i], misfits2[i], iters[i], diffs[i], true);

  /* 		int time, misfit, diff; 
      for (int j = 0; j < 4; j++) {  
        time 
        run_test(stdout, num, time, misfit, diff);
        times[i] += time;
        misfits[i] += misfit;
        diffs[i] += diff;
      } */

    }

    printf("\n------------ Results summary --------------------\n\n");
    printf("------------------------------------------------------------------------\n");
    printf("|  # |       Time | Iter |     Minimum |   Min(comp) |            Diff |\n");
    printf("------------------------------------------------------------------------\n");
    for( i = 1; i <= TEST_NUM; i++ ) {
      printf("| %2d | %10.6f | %4d | %11.7f | %11.7f | %1.13f |\n", i, times[i], misfits[i], misfits2[i], iters[i], diffs[i]);
    }
    printf("------------------------------------------------------------------------\n\n"); 
    
    

  }

  return(0);
}
