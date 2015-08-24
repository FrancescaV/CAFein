double signum(double);
double findMin(double, double );
void calcSteffenInterp(int,double *, double *,double **);
void evalSteffenInterp(int ,double *, double **, double ,double *);
void calcSteffenInterp_RiccatiElmnts(int, int, int, double **, double **, double **, double *, double **, double **);
void evalSteffenInterp_RiccatiElmnts(int, int, int, double *, double **, double, double *, gsl_matrix*, int);
