int func_Riccati_R(double csi, const double y[], double f[], void *params);
int jac_Riccati_R(double csi, const double y[], double *dfdy, double dfdt[], void *params);
int func_Riccati_Ronly(double csi, const double y[], double f[], void *params);
int jac_Riccati_Ronly(double csi, const double y[], double *dfdy, double dfdt[], void *params);
int func_Riccati_Vonly(double csi, const double y[], double f[], void *params);
int func_Riccati_Vonly_rkf45 (double csi, const double y[], double f[], void *params);
int jac_Riccati_Vonly(double csi, const double y[], double *dfdy, double dfdt[], void *params);
