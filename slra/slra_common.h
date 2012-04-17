         


void tmv_prod_new(gsl_matrix*, int, 
	      gsl_vector*, int, gsl_vector*);

/*int tls(gsl_matrix*, gsl_matrix*, gsl_matrix*);*/


void print_mat(const gsl_matrix*);
void print_mat_tr(const gsl_matrix*);
void print_arr(double*, int);

void gsl_matrix_vectorize(double*, gsl_matrix*);
void gsl_matrix_vec_inv(gsl_matrix*, double*);

void m_to_gsl_matrix(gsl_matrix* a_gsl, double* a_m);
void gsl_to_m_matrix(double* a_m, gsl_matrix* a_gsl); 

/* Convert double matrix to structure */
/*int slraMatrix2Struct( data_struct *s, double *s_matr, 
                       int q, int s_matr_cols );*/
void slraString2Method( const char *str_buf, opt_and_info *popt );
int slraString2Disp( const char *str_value );


