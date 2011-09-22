/* 
MEX function for calling the SLICOT function MB02GD followed by DPBTRS.
Solves a positive definite, block banded, Toeplitz system of equations. 

In Linux and Solaris execute: 

g77 -c MA02FD.f  MB02CU.f  MB02CV.f  MB02GD.f
ar -r LIBMB02GD.a MA02FD.o  MB02CU.o  MB02CV.o  MB02GD.o

Link to LAPACK and BLAS. Mex command (you will need to adapt the paths):

Solaris: mex mb02gd.c LIBMB02GD.a /software/lib/liblapack.a /software/lib/libblas.a 
Linux:   mex mb02gd.c LIBMB02GD.a /usr/lib/liblapack.a /usr/lib/libblas.a /usr/lib/gcc-lib/i386-redhat-linux/3.2/libg2c.a

Calling sequence: mb02gd(F,r), where F is the nonzero part of the first block row 
and r is the right-hand-side. 
*/

#include <stdio.h>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )

{
  double *f, *f_copy, *r, *rb, *dwork, *x;
  int    k, n, s, nl, fs, rs, ldwork, d, info, kd, zero = 0; 
  char   *chr = "R", *chn = "N", *chu = "U", msg[50];

  /* Get input data */
  f  = mxGetPr( prhs[0] ); 
  r  = mxGetPr( prhs[1] );

  /* Define constants */
  k  = mxGetM( prhs[0] );
  fs = mxGetN( prhs[0] );
  rs = mxGetM( prhs[1] );
  d  = mxGetN( prhs[1] );
  n  = rs / k;
  nl = fs / k - 1;
  kd = fs - 1;
  ldwork = 1000; /* * ( (nl + 1)*k*k + (3 + nl)*k ); */

  /* Create output and work arrays */
  plhs[0] = mxCreateDoubleMatrix(rs, d, mxREAL);
  x       = mxGetPr( plhs[0] );  
  rb      = mxCalloc(fs*rs, sizeof(double));
  dwork   = mxCalloc(ldwork, sizeof(double));
  f_copy  = mxCalloc(k*fs, sizeof(double));

  if ( x == NULL || rb == NULL || dwork == NULL || f_copy == NULL ) {
     mexErrMsgTxt("Not enough memory!\n");
  }

  /* Copy r to x and f to f_copy (to preserve f unchanged) */
  memcpy(x,r,rs*d*sizeof(double));
  memcpy(f_copy,f,k*fs*sizeof(double));

  /* Call MB02GD */ 
  mb02gd_(chr, chn, &k, &n, &nl, &zero, &n, f_copy, &k, rb, &fs, dwork, &ldwork, &info);
  
  if (info != 0) {
     sprintf(msg, "MB02GD failed with an error: info = %d.\n", info);
     mexErrMsgTxt(msg);
  }
  mxFree(dwork);
  mxFree(f_copy);

  /* Call DPBTRS */ 
  dpbtrs_(chu, &rs, &kd, &d, rb, &fs, x, &rs, &info);
  mxFree(rb);

  if (info != 0) {
     sprintf(msg, "DPBTRS failed with an error: info = %d.\n", info);
     mexErrMsgTxt(msg);
  }  
}

