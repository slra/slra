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
#include "slra.h"



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

void run_test( FILE * log, char * testname, double & time, double & fmin, double &fmin2, int & iter, double &diff, 
               char * method = "l", int use_slicot = 0, bool silent = false ) {
  gsl_matrix *xt = NULL, *x = NULL, *a = NULL, *b = NULL, *v = NULL, *perm = NULL;
  gsl_vector *p = NULL, * p2 = NULL;
  
  data_struct s; /* {1,2,{'T',10,1,'U',1,1}}; */
  opt_and_info opt;
  
  slraAssignDefOptValues(opt);
//  opt.maxiter = 500;
  opt.maxiter = 500;
  opt.disp = SLRA_OPT_DISP_ITER;
  slraString2Method(method, &opt);
  opt.use_slicot = use_slicot;
  
  
  int i, j, m = 9599, n = 12, d = 4, tmp, np = 9599;
  int x_given;
  char faname[MAX_FN_LEN], fbname[MAX_FN_LEN], fxname[MAX_FN_LEN], fpname[MAX_FN_LEN],
       fxtname[MAX_FN_LEN], fsname[MAX_FN_LEN], fxresname[MAX_FN_LEN],
       fpresname[MAX_FN_LEN];
  FILE *file;

  /* default file names */
/*  sprintf(faname,"a%s.txt",testname);
  sprintf(fbname,"b%s.txt",testname);*/
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
    fscanf(file, "%d %d %d %d %d %d", &m, &n, &d, &s.k, &s.q, &np); 
    if ((s.k <= 0) || (m % s.k != 0)) {
      fprintf(log, "Bad k: %d \n", s.k);
    }

    if ((s.q <= 0) || (s.q > 10)) {
	    fprintf(log, "Bad k: %d \n", s.k);
	    throw 1;
    }
    
    for (i = 0; i < s.q; i++)  {
      fscanf(file, "%d %d %d %d", &(s.a[i].blocks_in_row), &(s.a[i].nb), &(s.a[i].exact), &(s.a[i].toeplitz)); 
    }
    fclose(file);

    /* read the data a, b from a.txt and b.txt */
/*    if (((a = gsl_matrix_alloc(m, n)) == NULL) || !read_mat(a, faname, log)) {
      throw 1;
    }
    if (((b = gsl_matrix_alloc(m, d)) == NULL) || !read_mat(b, fbname, log)) {
      throw 1;
    }*/


    if (((x = gsl_matrix_calloc(n, d)) == NULL)) {
      throw 1;
    }

    if (((perm = gsl_matrix_calloc(n+d, n+d)) == NULL)) {
      throw 1;
    }

    
    x_given = read_mat(x, fxname, log);

    if (((p = gsl_vector_alloc(np)) == NULL) || !read_vec(p, fpname, log)) { 
      throw 1;
    }

    if (((p2 = gsl_vector_alloc(np)) == NULL)) { 
      throw 1;
    }
    
    gsl_vector_memcpy(p2,p);

    
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
    slra(p, &s, x, v, &opt, x_given, 1 /* Compute correction */, perm, 0);

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
    
//    fmin2 = opt.chol_time;
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
    if (perm != NULL) {
      gsl_matrix_free(perm);
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


#define TEST_NUM 7

int main(int argc, char *argv[])
{
  double times[TEST_NUM+1], misfits[TEST_NUM+1], misfits2[TEST_NUM+1],  diffs[TEST_NUM+1];
  int iters[TEST_NUM+1];
  char num[10];
  char * method = "l";
  int  i = -1;
  int use_slicot = 1;


  if (argc > 1) {
    method = argv[1];
  }

  if (argc > 3) {
    i = atoi(argv[3]);
  }

  if (argc > 2) {
    use_slicot = atoi(argv[2]);
  }


  if (i >= 1 && i <= TEST_NUM) {
    printf("\n------------------ Testing example %d  ------------------\n\n", i);
    sprintf(num, "%d", i);
    run_test(stdout, num, times[i], misfits[i], misfits2[i], iters[i], diffs[i], method, use_slicot);

    printf("\n------------ Results summary --------------------\n\n");
    printf("------------------------------------------------------------------------\n");
    printf("|  # |       Time | Iter |     Minimum | ||dp||^2    |        Diff (X) |\n");
    printf("------------------------------------------------------------------------\n");
    printf("| %2d | %10.6f | %4d | %11.7f | %11.7f | %15.10f |\n", i, times[i], misfits[i], misfits2[i], iters[i], diffs[i], method);
    printf("------------------------------------------------------------------------\n\n"); 
  } else {
    printf("\n------------------ Testing all examples  ------------------\n\n");

    for( i = 1; i <= TEST_NUM; i++ ) {
      sprintf(num, "%d", i);
      run_test(stdout, num, times[i], misfits[i], misfits2[i], iters[i], diffs[i], method, use_slicot);

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
    printf("|  # |       Time | Iter |     Minimum | ||dp||^2    |        Diff (X) |\n");
    printf("------------------------------------------------------------------------\n");
    for( i = 1; i <= TEST_NUM; i++ ) {
      printf("| %2d | %10.6f | %4d | %11.7f | %11.7f | %15.10f |\n", i, times[i], misfits[i], misfits2[i], iters[i], diffs[i], method);
    }
    printf("------------------------------------------------------------------------\n\n"); 
  }   
    


  return 0;
}
