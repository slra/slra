#include <memory.h>
#include "slra.h"

void a_kron_id( const gsl_matrix *A, size_t d,  gsl_matrix *IkronA ) {
  gsl_matrix_set_zero(IkronA);
  
  for (size_t  i = 0; i < A->size1; i++) {
    for (size_t  j = 0; j < A->size2; j++) {
      
      gsl_matrix subA = gsl_matrix_submatrix(IkronA, d * i, d * j, d, d).matrix;
      gsl_vector diag = gsl_matrix_diagonal(&subA).vector;
      gsl_vector_set_all(&diag, gsl_matrix_get(A, i, j));
    }
  }
}

NLSVarproPsiVecR::NLSVarproPsiVecR( VarproFunction &fun, const gsl_matrix *PsiT ) :
      NLSVarpro(fun) {
  if (PsiT == NULL) {
    myPsiTBig = gsl_matrix_alloc(myFun.getNrow() * myFun.getD(), myFun.getNrow() * myFun.getD());
    gsl_matrix_set_identity(myPsiTBig);
    myNEssVar = (myFun.getNrow() - myFun.getD()) * myFun.getD();
  } else {
    if (PsiT->size1 < PsiT->size2) {
      throw new Exception("Incorrect sizes of Psi matrix.\n");
    }

    if (PsiT->size1 == myFun.getNrow()) {
      myPsiTBig = gsl_matrix_alloc((PsiT->size1) * myFun.getD(), (PsiT->size2) * myFun.getD());
      a_kron_id(PsiT, myFun.getD(), myPsiTBig);
      myNEssVar = (PsiT->size2 - myFun.getD()) * myFun.getD();
    } else {
      myPsiTBig = gsl_matrix_alloc(PsiT->size1, PsiT->size2);
      gsl_matrix_memcpy(myPsiTBig, PsiT);
      myNEssVar = PsiT->size2;
    }
  }
        
  myTmpRVec = gsl_vector_alloc(myFun.getNrow() * myFun.getD());
  myTmpR = gsl_matrix_view_vector(myTmpRVec, myFun.getNrow(), myFun.getD()).matrix;
}

NLSVarproPsiVecR::~NLSVarproPsiVecR() {
  gsl_vector_free(myTmpRVec);
  gsl_matrix_free(myPsiTBig);
}

void NLSVarproPsiVecR::RTheta2x( const gsl_matrix *RTheta, gsl_vector *x )
{
  gsl_matrix_memcpy(&myTmpR, RTheta);
  
  gsl_matrix *At = gsl_matrix_alloc(myPsiTBig->size1, myPsiTBig->size2);
  gsl_matrix_memcpy(At, myPsiTBig);
  
  size_t one = 1, lwork = -1, info;
  double tmp;

  dgels_("T", &At->size2, &At->size1, &one, At->data,
         &At->size2, myTmpRVec->data, &At->size1, &tmp, &lwork, &info);
  double *work = new double[lwork = tmp];
  dgels_("T", &myPsiTBig->size2, &At->size1, &one, At->data,
         &At->size2, myTmpRVec->data, &At->size1, work, &lwork, &info);
  delete [] work;
  
  gsl_matrix_free(At);
  gsl_vector sol_x = gsl_vector_subvector(myTmpRVec, 0, x->size).vector;
  gsl_vector_memcpy(x, &sol_x);
}

void NLSVarproPsiVecR::x2RTheta( gsl_matrix *RTheta, const gsl_vector *x )
{
  gsl_blas_dgemv(CblasNoTrans, 1.0, myPsiTBig, x, 0.0, myTmpRVec);
  gsl_matrix_memcpy(RTheta, &myTmpR);
}


