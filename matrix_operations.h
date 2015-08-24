void matrix_mul_rows_by_cols (gsl_matrix *, gsl_matrix *, gsl_matrix *, int, gsl_vector *);
void matrix_sum (gsl_matrix *, gsl_matrix *, gsl_matrix *, int);
void matrix_diff (gsl_matrix *, gsl_matrix *, gsl_matrix *, int);
void extract_submatrices(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *,gsl_matrix *, int);
void mul_matrix_by_vec(gsl_matrix *, gsl_vector *, gsl_vector *, gsl_vector *, int);
double determinantNbyN(gsl_matrix *, int, gsl_permutation *, gsl_matrix *);
void complexDeterminant3by3(gsl_matrix *, gsl_matrix *, gsl_vector *);
void InverseMatrix(gsl_matrix *, gsl_matrix *, int, gsl_permutation *, gsl_matrix *);
void vector_sum (gsl_vector *, gsl_vector *, gsl_vector *, int );
