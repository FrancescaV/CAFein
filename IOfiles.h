int readInputFile(const char *, std::vector<std::vector<double> > &);
int countLinesInFile(const char *fileName);
int countElementsInInterval(std::ifstream &, int, double, double, int);
void fillInRiccatiElements(std::ifstream &, int, double, double, int, int, double *, double **, std::vector<std::vector<double> > &);
void OpenAndLabel_rkf45_state(const char *, int, std::ofstream &, int);
void OpenAndLabel_Rij_inOut(const char *, int, std::ofstream &, int, int);
void fillIn_rkf45_state(std::ifstream &, int, double, double, int, int, double **, std::vector<std::vector<double> > &);
void labelF_Riccati_xFit(std::ofstream &, int);
void writeF_Riccati_xFit(std::ofstream &, double, double, gsl_matrix *, gsl_matrix *, gsl_vector *, int, int);
void labelF_Riccati_inOut(std::ofstream &, int, int);
void writeF_Riccati_inOut(std::ofstream &, double, double, gsl_vector *, int, double, 
						  double, double, double, double, double, double, int);
void labelF_eigenfunctions(std::ofstream &, int, int ,int);
int writeF_eigenfunctions(std::ofstream &, double, double, double, double, double,  
						   double, double, gsl_vector *, gsl_vector *, gsl_matrix *, int, int, int, double, 
						   double, double, double, double, double, double, double, double);
void labelF_eigenfunctions_and_unperturbed_model(std::ofstream &, int, int, int);
int writeF_eigenfunctions_and_unperturbed_model(std::ofstream &, double, double, gsl_vector *, gsl_vector *, double, 
												double, double, double, double, double, double, double, double, double, 
												double, double, double, double, double, double, double, int, int, int, int, 
												double, double, double, double, double);
void labelF_globalProperties(std::ofstream &fileName);
void writeF_globalProperties(std::ofstream &fileName,  double,  double, double, double, double, double, double, double,  
							 double, double, double, double, double, double, double, double, double, double, double, double, double);
void labelF_integratorState_inOut(std::ofstream &, int);
void writeF_integratorState_inOut(std::ofstream &, gsl_vector *, int);
void errorMessageAndExit(const char *, const char *);
