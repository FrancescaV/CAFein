void create_Riccati_ABCD_new(double, double, double ,double ,double ,double ,double ,double, double, double, double, double, double,
							 double, double, double, double, double, double, double, double, 
							 gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, int, int, int, int, int);
void permute_Riccati_ABCD(double, double, double, double,double, double,double, double, double, double, double, double, double,
						  double, double, double, double, double, double, 
						  gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
						  gsl_matrix *,gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_permutation *, int, int, int, int, int);
void permute_Riccati_matrix_R(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix * ,gsl_matrix *, 
							  gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_vector *, 
							  gsl_permutation *, int);
void create_Riccati_eqs(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
						gsl_matrix *, int, gsl_matrix *, gsl_matrix *, gsl_vector *);
double calc_detRiccati_condition(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *,
								 gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
								 gsl_vector *, gsl_permutation *, int);
void back_to_original_Riccati_R(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *,
								gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
								gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_permutation *, int);
void find_permutation_with_minNorm(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
								   gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_vector *, 
								   gsl_vector *, gsl_vector *, gsl_permutation *, int, int &, int, int);
void fill_permutation_T(gsl_matrix *, gsl_vector *, int);
void permute_Riccati_UV(gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
						gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_vector *, 
						gsl_vector *, int);
void back_to_original_Riccati_UV(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
								 gsl_matrix *, gsl_vector *, gsl_vector *, gsl_vector *, 
								 gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *, 
								 gsl_vector *, gsl_vector *, gsl_permutation *, int);
void permute_Riccati_V(gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_matrix *, 
					   gsl_matrix *, gsl_matrix *, gsl_vector *, gsl_vector *, int);
void store_IC_ForIntegrationEigenfunctions(gsl_matrix *, double, double, gsl_vector *, gsl_vector *,
										   gsl_vector *, int &, int);
void update_csi_and_r_ij_forEndInterval(gsl_matrix *, double, gsl_vector *, int, int);
void store_Tstep_ForIntegrationEigenfunctions(gsl_matrix *, gsl_vector *, int, int &, int);
void store_Ttot_ForIntegrationEigenfunctions(gsl_matrix *, int, int, gsl_matrix *, int);
void create_permutedCD(double, double,double, double,double, double, double, double, double, double, double, 							
					   double, double, double, double, double, double, double, double, gsl_matrix *, gsl_matrix *, int, int);
void complete_y_Riccati_usingSimmetryR(double *, int);
void complete_matR_usingSymmetry(gsl_matrix *, int);
