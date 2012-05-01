         
/*
 * tmv_prod_new: block-Toeplitz banded matrix p =  T * v
 * T - storage for [t_s-1' ... t_1' t_0 t_1 ... t_s-1].
 * m = number of block rows / columns
 */ 
void tmv_prod_vector( gsl_vector *T, int s, gsl_vector* v, int m, 
         gsl_vector* p );
void tmv_prod_new( gsl_matrix *T, int s,  gsl_vector *v, int m, gsl_vector *p, 
         double beta = 0.0 );
         
void copyLowerTrg( gsl_matrix * dest, const gsl_matrix *src  );         
void shiftLowerTrg( gsl_matrix * dest, const gsl_matrix *src  );         

void print_mat(const gsl_matrix*);
void print_mat_tr(const gsl_matrix*);
void print_arr(const double*, int);

void print_vec(const gsl_vector*);

void gsl_matrix_vectorize(double*, gsl_matrix*);
void gsl_matrix_vec_inv(gsl_matrix*, double*);

void m_to_gsl_matrix(gsl_matrix* a_gsl, double* a_m);
void gsl_to_m_matrix(double* a_m, gsl_matrix* a_gsl); 

/* Convert double matrix to structure */
/*int slraMatrix2Struct( data_struct *s, double *s_matr, 
                       int q, int s_matr_cols );*/
void String2Method( const char *str_buf, opt_and_info *popt );
int String2Disp( const char *str_value );

Structure *createMosaicStructure( gsl_vector * ml,  gsl_vector *nk, 
               gsl_vector * wk, int np_comp );
               
int compute_np( gsl_vector* ml, gsl_vector *nk );               


const gsl_vector *vecChkNIL( const gsl_vector &vec );
gsl_vector *vecChkNIL( gsl_vector &vec  );
gsl_matrix *matChkNIL( gsl_matrix &mat_vw );

void tolowerstr( char * str );
