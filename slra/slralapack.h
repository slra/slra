 
#ifdef __cplusplus
extern "C" {
#endif


/* SLICOT and LAPACK functions */

#ifdef USE_SLICOT

void mb02gd_(const char *typet, const char *triu, const size_t *k, const size_t *n, const size_t *nl, 
             const size_t *p, const size_t *s, double *t, const size_t* ldt, 
             double *rb, const size_t *ldrb, double *dwork, const size_t *ldwork, size_t *info);
#endif             
             
void dtbtrs_(const char* uplo, const char* trans, const char* diag, const size_t* n, const  size_t* kd, const size_t* nrhs, 
             const double* ab, const size_t* ldab, const double* b, const size_t* ldb, size_t* info); 

void dpbtrs_(const char* uplo, const size_t* n, const size_t* kd, const size_t* nrhs, 
             const double* ab, const size_t* ldab, const double* b, const size_t* ldb,  size_t* info); 

void dpbtrf_(const char* uplo, const size_t* n, const size_t* kd, double* ab, const size_t* ldab, size_t* info); 

void dgelqf_(const size_t *m, const size_t *n, double *a, const  size_t *lda, 
             double *tau, double *work, const size_t *ldwork, size_t *info);
             
              
void dormlq_(const char *side, const char *trans, const size_t *m, const size_t *n, const size_t *k,
             double *a, const size_t *lda, double *tau, double *c, const size_t *ldc, 
             double *work, const size_t *lwork, size_t *info);
void dtrsm_(const char* side, const char *uplo, const char *transa, const char *diag,
            const size_t *m, const size_t *n, const double *alpha, const double *a, const size_t *lda,
            double *b, const size_t *ldb);


void dgesv_(const size_t* n, const size_t* nrhs, double* a, const size_t* lda, 
            const size_t* ipiv, double* b, const size_t* ldb, size_t* info);
            
void dgesvd_(const char* jobu, const char* jobvt, const size_t* m, const size_t* n, double* a, const size_t* lda, double* s, const double* u, const size_t* ldu, const double* vt, const size_t* ldvt, double* work, const size_t* lwork, size_t * info);              


#ifdef __cplusplus
}
#endif
