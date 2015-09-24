#include "slra.h"

PhiStructure::PhiStructure( const gsl_matrix *PhiT, Structure *S  ) :
  myPStruct(S) {
  if (PhiT->size1 < PhiT->size2) {
    throw new Exception("PhiStructure::PhiStructure");
  }
  myPhiT = gsl_matrix_alloc(PhiT->size1, PhiT->size2);
  gsl_matrix_memcpy(myPhiT, PhiT);
  myTempStMat = gsl_matrix_alloc(myPStruct->getN(), myPStruct->getM());
}

PhiStructure::~PhiStructure()  {
  if (myPStruct != NULL) {
    delete myPStruct;
  }
  gsl_matrix_free(myPhiT);
  gsl_matrix_free(myTempStMat);
}

void PhiStructure::fillMatrixFromP( gsl_matrix* c, const gsl_vector* p )
{
  myPStruct->fillMatrixFromP(myTempStMat, p);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myTempStMat, myPhiT, 0, c);
}

gsl_matrix *PhiStructure::createPhiTRt( const gsl_matrix *Rt ) const {
  gsl_matrix *res = gsl_matrix_alloc(myPhiT->size1, Rt->size2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, myPhiT, Rt, 0, res);
  return res;
}

void PhiStructure::multByGtUnweighted( gsl_vector* p, const gsl_matrix *Rt,
        const gsl_vector *y, double alpha, double beta, bool skipFixedBlocks ) {
  gsl_matrix *PhiTRt = createPhiTRt(Rt);
  myPStruct->multByGtUnweighted(p, PhiTRt, y, alpha, beta, skipFixedBlocks);
  gsl_matrix_free(PhiTRt);
}

PhiStructure::PhiCholesky::PhiCholesky( const PhiStructure *s, size_t d ) :
   myStruct(s), myParent(s->myPStruct->createCholesky(d)) {
}

PhiStructure::PhiCholesky::~PhiCholesky() {
  if (myParent != NULL) {
    delete myParent;
  }
}

void PhiStructure::PhiCholesky::calcGammaCholesky( const gsl_matrix *Rt, double reg ) {
  gsl_matrix *PhiTRt = myStruct->createPhiTRt(Rt);
  myParent->calcGammaCholesky(PhiTRt, reg);
  gsl_matrix_free(PhiTRt);
}

void PhiStructure::PhiCholesky::multInvCholeskyVector( gsl_vector * y_r, long trans ) {
  myParent->multInvCholeskyVector(y_r, trans);
}

void PhiStructure::PhiCholesky::multInvGammaVector( gsl_vector * y_r ) {
  myParent->multInvGammaVector(y_r);
}

PhiStructure::PhiDGamma::PhiDGamma( const PhiStructure *s, size_t d ) :
    myStruct(s), myParent(s->myPStruct->createDGamma(d)) {
  myPhi = gsl_matrix_alloc(myStruct->myPhiT->size2, myStruct->myPhiT->size1);
  gsl_matrix_transpose_memcpy(myPhi, myStruct->myPhiT);
}

PhiStructure::PhiDGamma::~PhiDGamma() {
  if (myParent != NULL) {
    delete myParent;
  }
  if (myPhi != NULL) {
    gsl_matrix_free(myPhi);
  }
}

void PhiStructure::PhiDGamma::calcYtDgammaY( gsl_matrix *At,
                           const gsl_matrix *Rt, const gsl_matrix *Yt ) {
  gsl_matrix *PhiTRt = myStruct->createPhiTRt(Rt);
  gsl_matrix *tempAt = gsl_matrix_alloc(PhiTRt->size1, PhiTRt->size2);
  myParent->calcYtDgammaY(tempAt, PhiTRt, Yt);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, myStruct->myPhiT, tempAt, 0, At);
  gsl_matrix_free(PhiTRt);
  gsl_matrix_free(tempAt);
}

void PhiStructure::PhiDGamma::calcDijGammaYr( gsl_vector *z,
         const gsl_matrix *Rt, size_t j_1, size_t i_1, const gsl_vector *y,
         const gsl_matrix *Phi ) {
  if (Phi != NULL) {
    throw new Exception("Does not support nested Phi multiplication...");
  }
  gsl_matrix *PhiTRt = myStruct->createPhiTRt(Rt);
  myParent->calcDijGammaYr(z, PhiTRt, j_1, i_1, y, myPhi);
  gsl_matrix_free(PhiTRt);
}



