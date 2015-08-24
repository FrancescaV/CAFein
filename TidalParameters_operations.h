void calculateTidalCouplingCoeff(int, double *, double **, double *, int , double *, double *, 
								 double *, double *, double *, double, double, double, double *);
double calculate_LegendrePol_Plm0(int, int);
double calculateFourierCoeffs_Clmk(int, int, int, double, double, double);
double calculate_Glmk_2(int, int, int, double, double, double);
double calculate_Glmk_3(int, int, int, double, double, double);
void calculate_Flmk(int, int, int, double, double, double, double, double, double, gsl_vector *, gsl_vector *);
void anomalyGrid(double, gsl_vector *, gsl_vector *);
double calculateCoeffs_Klmk(int, int);
double integrate_with_extended_trapezium_rule (gsl_vector* , gsl_vector* , int); 
double calculate_dAdT_tides(int, int, int, double, double, double, double, double, double, double, int, double, double);
double calculate_dAdT_GR(double, double, double, double);
double calculate_dOmegadT_tides(int, int, int, double, double, double, double, double, double, double, int, double, double, double);
double calculate_dedT_tides(int, int, int, double, double, double, double, double, double, double, int, double, double);


