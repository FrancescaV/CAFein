#include <vector>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include "SteffenInterp.h"
#include "TidalParameters_operations.h"
#include "IOfiles.h"
#include "numericalIntegration.h"
#include "dma.h"

using namespace std;

const double PI = 2.0*acos(0.0);

const int ngrid = 1e6;

/*************************************************************************************
 ************ 
 ************	The routines described below handle the tidal parameters
 ************	Q(tidal overlap integral), c_lmk (Fourier coefficients)...
 ************ 
 **************************************************************************************/

/**************************************************************************************
 **** 01/25/2013 I checked the accuracy of Clmk, G^(2)_lmk, and G^(3)_lmk
 **** by re-calculating them using the new integrator from
 **** Gill and Miller and the numbers change only at the ~13th decimal place.
 **************************************************************************************/



/*************************************************************************************
 ****** This routine is used to calculated the tidal coupling coefficient
 ****** following Eq. 5.22 in Bart's thesis.
 ****** 
 ****** To calculate the integral the trapezium rule is used:
 ****** int from A to B of f(x)dx = (B-A)*[f(B)+F(A)]/2
 ****** 
 ****** Note that: rhoRlP1 = rho*r^(l+1) 
 ****** 
 ****** Input: 
 ****** - meshModel_func = number of mesh points in the stellar model
 ****** - radius_model_func = vector containing the model's radius
 ****** - fitRhoCoeffs_func = Steffen interpolation's coefficients for the model's density rho
 ****** - rhoFuncs_func = vector containing the result of Steffen interpolant
 ****** - meshEigenf_func = number of mesh points in the eigenfunction
 ****** - csi_Tr_func = vector containing the radial component of the displacement
 ****** - csi_Th_func = vector containing the orthogonal component of the displacement
 ****** - radiusEigenf_func = vector containing the eigenfunction's radius
 ****** - deg_l_func = l-order of the mode
 ****** - w_r_func = real eigenfrequency
 ****** - w_i_func = imaginary eigenfrequency
 ******
 ****** Output
 ****** - Q_tidalCoupling_func = tidal coupling coefficient (real and imaginary components). 
 ******  
 ****** Another way of calculating Q_tidCoup is given by:
 ****** Press & Teukolsky 1977 APJ 213: 183-192 Eq. (38)
 ****** To chose this option uncomment the part below.
 **************************************************************************************/
void calculateTidalCouplingCoeff(int meshModel_func, double *radius_model_func, double **fitRhoCoeffs_func, double *rhoFuncs_func, 
								 int meshEigenf_func, double *csi_Tr_real_func, double *csi_Tr_imag_func, 
								 double *csi_Th_real_func, double *csi_Th_imag_func, double *radiusEigenf_func, double deg_l_func, 
								 double w_r_func, double w_i_func, double *Q_tidalCoupling_func)

{
	double Ns_rad = 0.0, Ns_orth = 0.0, Ns = 0.0, rhoR2_A = 0.0, rhoR2_B = 0.0, fA = 0.0, fB = 0.0, delta_r = 0.0, 
	Qrad_real = 0.0, Qrad_imag = 0.0, Qorth_real = 0.0, Qorth_imag = 0.0, rhoRlP1_A = 0.0, rhoRlP1_B = 0.0, rA = 0.0, rB = 0.0, 
	Q_real = 0.0, Q_imag = 0.0;
	
	cout << "meshEigenf_func = " << meshEigenf_func << endl;
	for (int i=1; i< meshEigenf_func; i++)
	{				
		
		rA = radiusEigenf_func[i-1];
		rB = radiusEigenf_func[i];		
		delta_r = rB - rA;
		
		evalSteffenInterp(meshModel_func,radius_model_func,fitRhoCoeffs_func,rA, rhoFuncs_func);			
		rhoR2_A = rhoFuncs_func[0]*rA*rA; //rho r^2
		rhoRlP1_A = rhoFuncs_func[0]*pow(rA, deg_l_func + 1.0); //rho r^(l+1)
				
		evalSteffenInterp(meshModel_func,radius_model_func,fitRhoCoeffs_func,rB, rhoFuncs_func);			
		rhoR2_B = rhoFuncs_func[0]*rB*rB;
		rhoRlP1_B = rhoFuncs_func[0]*pow(rB, deg_l_func + 1.0);
		
		/*for adiabatic equations the imaginary csi_T is null --> it falls back to the adiabatic one*/
		fA = (csi_Tr_real_func[i-1]*csi_Tr_real_func[i-1]+csi_Tr_imag_func[i-1]*csi_Tr_imag_func[i-1])*rhoR2_A;
		fB = (csi_Tr_real_func[i]*csi_Tr_real_func[i]+csi_Tr_imag_func[i]*csi_Tr_imag_func[i])*rhoR2_B;		
		Ns_rad = Ns_rad + delta_r*(fB + fA)/2.0;
		
		/*for adiabatic equations the imaginary csi_T is null --> it falls back to the adiabatic one*/
		fA = (csi_Th_real_func[i-1]*csi_Th_real_func[i-1]+csi_Th_imag_func[i-1]*csi_Th_imag_func[i-1])*rhoR2_A;
		fB = (csi_Th_real_func[i]*csi_Th_real_func[i]+csi_Th_imag_func[i]*csi_Th_imag_func[i])*rhoR2_B;
		Ns_orth = Ns_orth + delta_r*(fB + fA)/2.0;
				
		fA = csi_Tr_real_func[i-1]*rhoRlP1_A;
		fB = csi_Tr_real_func[i]*rhoRlP1_B;		
		Qrad_real = Qrad_real + delta_r*(fB + fA)/2.0;

		/*for adiabatic equations the imaginary csi_Tr is null Qrad_imag is null*/
		fA = csi_Tr_imag_func[i-1]*rhoRlP1_A;
		fB = csi_Tr_imag_func[i]*rhoRlP1_B;		
		Qrad_imag = Qrad_imag + delta_r*(fB + fA)/2.0;
				
		fA = csi_Th_real_func[i-1]*rhoRlP1_A;
		fB = csi_Th_real_func[i]*rhoRlP1_B;		
		Qorth_real = Qorth_real + delta_r*(fB + fA)/2.0;
		
		/*for adiabatic equations the imaginary csi_Th is null Qorth_imag is null*/
		fA = csi_Th_imag_func[i-1]*rhoRlP1_A;
		fB = csi_Th_imag_func[i]*rhoRlP1_B;		
		Qorth_imag = Qorth_imag + delta_r*(fB + fA)/2.0;
				
	}//for (i=0; i< meshEigenf_func; i++)
	Ns_orth = deg_l_func*(deg_l_func+1.0)*Ns_orth; 
	Ns = Ns_rad+Ns_orth;
	
	Qrad_real = deg_l_func*Qrad_real;
	Qorth_real = deg_l_func*(deg_l_func + 1.0)*Qorth_real;
	Q_real = Qrad_real + Qorth_real;
	
	/*for adiabatic equations this is null*/
	Qrad_imag = deg_l_func*Qrad_imag;
	Qorth_imag = deg_l_func*(deg_l_func + 1.0)*Qorth_imag;
	Q_imag = Qrad_imag + Qorth_imag;

	cout << "Ns in function calculateTidalCouplingCoeff is = " << Ns << endl;
	
	/*real tidal coupling:*/
	Q_tidalCoupling_func[0] = (1.0/Ns)*(1.0/((w_r_func*w_r_func + w_i_func*w_i_func)*(w_r_func*w_r_func + w_i_func*w_i_func)))*((w_r_func*w_r_func - w_i_func*w_i_func)*Q_real + 2.0*w_r_func*w_i_func*Q_imag);

	/*imaginary tidal coupling:*/
	Q_tidalCoupling_func[1] = (1.0/Ns)*(1.0/((w_r_func*w_r_func + w_i_func*w_i_func)*(w_r_func*w_r_func + w_i_func*w_i_func)))*((w_r_func*w_r_func - w_i_func*w_i_func)*Q_imag - 2.0*w_r_func*w_i_func*Q_real);
	
	//	/* Following Press & Teukolsky 1977 APJ 213: 183-192 Eq. (38)*/
	//	double Ns_rad = 0.0, Ns_orth = 0.0, Ns = 0.0, rhoR2_A = 0.0, rhoR2_B = 0.0, fA = 0.0, fB = 0.0, delta_r = 0.0, 
	//	Qrad = 0.0, Qorth = 0.0, Q_tidCoup = 0.0, rhoRlP1_A = 0.0, rhoRlP1_B = 0.0, rA = 0.0, rB = 0.0 ;
	//	
	//	for (int i=1; i< meshEigenf_func; i++)
	//	{						
	//		rA = radiusEigenf_func[i-1];
	//		rB = radiusEigenf_func[i];
	//		
	//		delta_r = rB - rA;
	//		
	//		evalSteffenInterp(meshModel_func,radius_model_func,fitRhoCoeffs_func,rA, rhoFuncs_func);			
	//		rhoRlP1_A = rhoFuncs_func[0]*pow(rA, deg_l_func + 1.0);
	//		
	//		evalSteffenInterp(meshModel_func,radius_model_func,fitRhoCoeffs_func,rB, rhoFuncs_func);			
	//		rhoRlP1_B = rhoFuncs_func[0]*pow(rB, deg_l_func + 1.0);
	//		
	//		fA = csi_Tr_real_func[i-1]*rhoRlP1_A;
	//		fB = csi_Tr_real_func[i]*rhoRlP1_B;
	//		
	//		Qrad = Qrad + delta_r*(fB + fA)/2.0;
	//		
	//		fA = csi_Th_real_func[i-1]*rhoRlP1_A;
	//		fB = csi_Th_real_func[i]*rhoRlP1_B;
	//		
	//		Qorth = Qorth + delta_r*(fB + fA)/2.0;
	//		
	//	}//for (i=0; i< meshEigenf_func; i++)
	//	
	//	Qrad = deg_l_func*Qrad;
	//	Qorth = deg_l_func*(deg_l_func + 1.0)*Qorth;
	//	Q_tidCoup = (Qrad + Qorth);
}


/*************************************************************************************
 ****** This routine is used to calculated the Legendre Polynimials P^|m|_l(0)
 ****** following Eq. 1.37 in Bart's thesis
 ****** 
 ****** Note that: 
 ****** 
 ****** - Gamma_LpMp1 = Gamma(l + |m| + 1)
 ****** - Gamma_LmMp2 = Gamma(l - |m| + 2)
 ****** 
 **************************************************************************************/

double calculate_LegendrePol_Plm0(int deg_l_int_func, int deg_m_int_func)
{
			
	double Plm0 = 0.0, Gamma_LpMp1 = 0.0, Gamma_LmMp2 = 0.0, absM = 0.0,
	deg_l = static_cast<double>(deg_l_int_func),
	deg_m = static_cast<double>(deg_m_int_func);
	
	
	
	if (2*((deg_l_int_func+deg_m_int_func)/2) != deg_l_int_func+deg_m_int_func) 
	{
		Plm0 = 0.0;
		return Plm0;
	}
	
	
	absM = fabs(deg_m);
	
	Gamma_LpMp1 = gsl_sf_gamma((deg_l + absM + 1.0)/2.0);
	
	Gamma_LmMp2 = gsl_sf_gamma((deg_l - absM + 2.0)/2.0);
	
	
	Plm0 = (pow(2.0, absM)/sqrt(PI))*(Gamma_LpMp1/Gamma_LmMp2)*cos((deg_l+absM)*PI/2.0);
	
	
	return Plm0;
}//calculate_LegendrePol_Plm0(double deg_l_func, double deg_m_func)

/*************************************************************************************
 ****** This routine is used to calculated the Fourier coefficients Clmk 
 ****** for an eccentric orbit and it is taken from Bart Willems' code Triton
 ****** C_lmk are calculated using Eq. 4 in Willems et al.(2010) 2010ApJ...713..239W
 ****** 
 ****** Note that: 
 ****** 
 ****** - Gamma_LmM = Gamma(l - |m|)
 ****** - Gamma_LpM = Gamma(l + |m|)
 ****** - Plm0 = associated Legendre polynomial of the first kind P^|m|_l(0)
 ****** 
 **************************************************************************************/
double calculateFourierCoeffs_Clmk(int deg_l_int_func, int deg_m_int_func, int deg_k_int_func, 
											 double ecc_func, double R1_noDim_func, double a_noDim_func)
{
	int i = 0;
	/****************************************************************
	 * Special cases: clmk = 0 when l+|m| is odd (because the associated
	 * Legendre function of the first kind evaluated at x=0 is equal to
	 * zero in this case), when k=0 and |m| > l-1 (this follows from																									   
	 * the binomial theorem), and when e=0 and k /= -m.
	 ****************************************************************/	
	double c_lmk = 0.0, Plm0 = 0.0, trueAnVal = 0.0, meanAnVal = 0.0,
	integralVal = 0.0, integrandVal = 0.0, absM = 0.0, Gamma_LmM = 0.0, Gamma_LpM = 0.0,
	deg_l = static_cast<double>(deg_l_int_func), 
	deg_m = static_cast<double>(deg_m_int_func), 
	deg_k = static_cast<double>(deg_k_int_func);
	
	absM = fabs(deg_m);
	
	if ((fabs(ecc_func - 0.0) <= 1.0e-16) && (deg_k_int_func != -deg_m_int_func))
	{
		c_lmk = 0.0;
		return c_lmk;
	}
	
	if (2*((deg_l_int_func+deg_m_int_func)/2) != deg_l_int_func+deg_m_int_func)
	{
		c_lmk = 0.0;
		return c_lmk;
	}
	
	if ((deg_k_int_func == 0) && (fabs(deg_m_int_func) > deg_l_int_func-1))
	{
		c_lmk = 0.0;
		return c_lmk;
	}
	
	
	/****************************************************************
	 **** Get arrays with grids of true and mean anomalies describing half
	 **** an orbit (i.e. true and mean anomaly from 0 to pi)
	 ****************************************************************/	
	gsl_vector *trueAn		= gsl_vector_alloc(ngrid);
	gsl_vector *meanAn		= gsl_vector_alloc(ngrid);
	gsl_vector *integrand	= gsl_vector_alloc(ngrid);
		
	anomalyGrid(ecc_func, trueAn, meanAn);
	
	/***************************************************************
	 **** Calculate the integrand of clmk using the trapezium rule
	 ***************************************************************/	
	for(i = 0; i < ngrid; i++)
	{
		trueAnVal = gsl_vector_get(trueAn, i);
		meanAnVal = gsl_vector_get(meanAn, i);
		integrandVal = cos(deg_m*trueAnVal+deg_k*meanAnVal)*
		pow(1.0+ecc_func*cos(trueAnVal), (deg_l-1.0));
		gsl_vector_set(integrand, i, integrandVal);
		
	}//for(i = 1; i < ngrid; i++)
	
	integralVal = integrate_with_extended_trapezium_rule (trueAn, integrand, ngrid);
	
	/****************
	 * free vectors
	 ****************/	
	gsl_vector_free(integrand);
	gsl_vector_free(meanAn);
	gsl_vector_free(trueAn);
	
	
	Plm0 = calculate_LegendrePol_Plm0(deg_l_int_func, deg_m_int_func);
	Gamma_LmM = gsl_sf_gamma(deg_l - absM + 1.0);
	Gamma_LpM = gsl_sf_gamma(deg_l + absM + 1.0);
	
	c_lmk = integralVal*(Gamma_LmM/Gamma_LpM)*(Plm0/PI)*pow(R1_noDim_func/a_noDim_func, deg_l - 2.0)*
	(1.0/pow(1.0 - ecc_func*ecc_func, deg_l - 0.5));
		
	return c_lmk;
}


/*************************************************************************************
 ****** This routine is used to calculated the Fourier coefficients Glmk_2
 ****** for an eccentric orbit and it is taken from Bart Willems' code Triton
 ****** Glmk_2 are calculated using Eq. 56 in Willems et al.(2010) 2010ApJ...713..239W
 **************************************************************************************/
double calculate_Glmk_2(int deg_l_int_func, int deg_m_int_func, int deg_k_int_func,
								  double ecc_func, double R1_noDim_func, double a_noDim_func)
{
	int i = 0;
	double Clmk = 0.0, Plm0 = 0.0, Glmk_2 = 0.0, integrandVal = 0.0, 
	integralValLeft = 0.0, integralValRight = 0.0, trueAnVal = 0.0, meanAnVal = 0.0, 
	deg_l = static_cast<double>(deg_l_int_func),
	deg_m = static_cast<double>(deg_m_int_func),
	deg_k = static_cast<double>(deg_k_int_func);
	
	
	Clmk = calculateFourierCoeffs_Clmk(deg_l_int_func, deg_m_int_func, deg_k_int_func, 
												 ecc_func, R1_noDim_func, a_noDim_func);

	/************************************************************************************* 
	 * Special cases: G2lmk = 0 when l+|m| is odd, when clmk = 0,
	 * and when m=k=0
	 **************************************************************************************/
	if (2*((deg_l_int_func+deg_m_int_func)/2) != deg_l_int_func+deg_m_int_func)
	{
		Glmk_2 = 0.0;
		return Glmk_2;
	}
	else if ((deg_m_int_func == 0) && (deg_k_int_func == 0))
	{
		Glmk_2 = 0.0;
		return Glmk_2;
	}
	else if (fabs(Clmk - 0.0) <= 1.0e-16)
	{
		Glmk_2 = 0.0;
		return Glmk_2;
	}
	
	/****************************************************************
	 * Get arrays with grids of true and mean anomalies describing half
	 * an orbit (i.e. true and mean anomaly from 0 to pi)
	 ****************************************************************/	
	gsl_vector *trueAn		= gsl_vector_alloc(ngrid);
	gsl_vector *meanAn		= gsl_vector_alloc(ngrid);
	gsl_vector *integrand	= gsl_vector_alloc(ngrid);
	
	anomalyGrid(ecc_func, trueAn, meanAn);
	
	/***************************************************************
	 * Calculate the integrand of clmk using the trapezium rule
	 ***************************************************************/	
	for(i = 0; i < ngrid; i++)
	{
		trueAnVal = gsl_vector_get(trueAn, i);
		meanAnVal = gsl_vector_get(meanAn, i);
		
		integrandVal = sin(deg_m*trueAnVal+deg_k*meanAnVal)*
		sin(trueAnVal)*pow(1.0+ecc_func*cos(trueAnVal), deg_l);		
		
		gsl_vector_set(integrand, i, integrandVal);
		
	}//for(i = 1; i < ngrid; i++)
	
	integralValLeft = integrate_with_extended_trapezium_rule (trueAn, integrand, ngrid);
	
	integralValLeft = (deg_l+1.0)*ecc_func*integralValLeft;
	
	for(i = 0; i < ngrid; i++)
	{
		trueAnVal = gsl_vector_get(trueAn, i);
		meanAnVal = gsl_vector_get(meanAn, i);
		
		integrandVal = cos(deg_m*trueAnVal+deg_k*meanAnVal)
		*pow(1.0+ecc_func*cos(trueAnVal), deg_l+1.0);	
		
		gsl_vector_set(integrand, i, integrandVal);
		
	}//for(i = 1; i < ngrid; i++)
	
	//tu multiply to left integral 
	integralValRight = integrate_with_extended_trapezium_rule (trueAn, integrand, ngrid);
	
	
	integralValRight = -deg_m*integralValRight;
	
	
	/*****************
	 * free vectors
	 *****************/	
	gsl_vector_free(integrand);
	gsl_vector_free(meanAn);
	gsl_vector_free(trueAn);
	
	
	Plm0 = calculate_LegendrePol_Plm0(deg_l_int_func, deg_m_int_func);	
	
	
	Glmk_2 = (2.0/pow(1.0-ecc_func*ecc_func, deg_l + 1.0))*Clmk*(Plm0/PI)*
	(integralValLeft + integralValRight);
	
	return Glmk_2;
}//double calculate_Glmk_2


/*************************************************************************************
 ****** This routine is used to calculated the Fourier coefficients Glmk_3
 ****** for an eccentric orbit and it is taken from Bart Willems' code Triton
 ****** Glmk_3 are calculated using Eq. 57 in Willems et al.(2010) 2010ApJ...713..239W
 **************************************************************************************/
double calculate_Glmk_3(int deg_l_int_func, int deg_m_int_func, int deg_k_int_func, 
								  double ecc_func, double R1_noDim_func, double a_noDim_func)
{
	int i = 0;
	double Clmk = 0.0, Plm0 = 0.0, Glmk_3 = 0.0, integrandVal = 0.0, 
	integralValLeft = 0.0, integralValRight = 0.0, trueAnVal = 0.0, meanAnVal = 0.0, 
	deg_l = static_cast<double>(deg_l_int_func),
	deg_m = static_cast<double>(deg_m_int_func),
	deg_k = static_cast<double>(deg_k_int_func);
	
	
	Clmk = calculateFourierCoeffs_Clmk(deg_l_int_func, deg_m_int_func, deg_k_int_func, 
												 ecc_func, R1_noDim_func, a_noDim_func);
	
	/************************************************************************************* 
	 * Special cases: G3lmk = 0 when l+|m| is odd, when clmk = 0,
	 * when e=0, and when m=k=0
	 **************************************************************************************/
	if (2*((deg_l_int_func+deg_m_int_func)/2) != deg_l_int_func+deg_m_int_func)
	{
		Glmk_3 = 0.0;
		return Glmk_3;
	}
	else if ((deg_m_int_func == 0) && (deg_k_int_func == 0))
	{
		Glmk_3 = 0.0;
		return Glmk_3;
	}
	else if (Clmk == 0.0)
	{
		Glmk_3 = 0.0;
		return Glmk_3;
	}
	else if (ecc_func == 0.0) 
	{
		Glmk_3 = 0.0;
		return Glmk_3;
	}
	/****************************************************************
	 * Get arrays with grids of true and mean anomalies describing half
	 * an orbit (i.e. true and mean anomaly from 0 to pi)
	 ****************************************************************/	
	gsl_vector *trueAn		= gsl_vector_alloc(ngrid);
	gsl_vector *meanAn		= gsl_vector_alloc(ngrid);
	gsl_vector *integrand	= gsl_vector_alloc(ngrid);
	
	anomalyGrid(ecc_func, trueAn, meanAn);
		
	/***************************************************************
	 * Calculate the integrand of clmk using the trapezium rule
	 ***************************************************************/	
	for(i = 0; i < ngrid; i++)
	{
		trueAnVal = gsl_vector_get(trueAn, i);
		meanAnVal = gsl_vector_get(meanAn, i);
		
		integrandVal = sin(deg_m*trueAnVal+deg_k*meanAnVal)*
		sin(trueAnVal)*pow(1.0+ecc_func*cos(trueAnVal), deg_l);		
		
		gsl_vector_set(integrand, i, integrandVal);
		
	}//for(i = 1; i < ngrid; i++)
	
	integralValLeft = integrate_with_extended_trapezium_rule (trueAn, integrand, ngrid);
	
	integralValLeft = (deg_l+1.0)*ecc_func*integralValLeft;
	
	for(i = 0; i < ngrid; i++)
	{
		trueAnVal = gsl_vector_get(trueAn, i);
		meanAnVal = gsl_vector_get(meanAn, i);
		
		integrandVal = cos(deg_m*trueAnVal+deg_k*meanAnVal)
		*pow(1.0+ecc_func*cos(trueAnVal), deg_l-1.0)*
		(pow(1.0+ecc_func*cos(trueAnVal), 2.0) - (1.0 - ecc_func*ecc_func));	
		
		gsl_vector_set(integrand, i, integrandVal);
		
	}//for(i = 1; i < ngrid; i++)
	
	//tu multiply to left integral 
	integralValRight = integrate_with_extended_trapezium_rule (trueAn, integrand, ngrid);
	
	integralValRight = -deg_m*integralValRight;
		
	/******************
	 * free vectors
	 ******************/
	gsl_vector_free(integrand);
	gsl_vector_free(meanAn);
	gsl_vector_free(trueAn);
		
	Plm0 = calculate_LegendrePol_Plm0(deg_l_int_func, deg_m_int_func);		
	
	Glmk_3 = (1.0/(ecc_func*pow(1.0-ecc_func*ecc_func, deg_l)))*Clmk*(Plm0/PI)*
	(integralValLeft + integralValRight);
	
	return Glmk_3;
}//double calculate_Glmk_3


/*************************************************************************************
 ****** This routine is used to calculated the Fourier coefficients Flmk
 ****** and it is taken from Bart Willems' code Triton
 ****** Flmk are calculated using Eq. 50 in Willems et al.(2010) 2010ApJ...713..239W
 **************************************************************************************/
void calculate_Flmk(int deg_l_int_func, int deg_m_int_func, int deg_k_int_func, 
					double ecc_func, double R1_noDim_func, double a_noDim_func, 						
					double epsilon_T_noDim_func, double g_func, double csiRel_func, 
					gsl_vector *Flmk_realImag_func, gsl_vector *Psi_T_realImag_func)
{

	double Clmk = 0.0, FlmkReal = 0.0, FlmkImag = 0.0, y3Real = 0.0, y3Imag = 0.0;
		
	Clmk = calculateFourierCoeffs_Clmk(deg_l_int_func, deg_m_int_func, deg_k_int_func, ecc_func, R1_noDim_func, a_noDim_func);

	y3Real = gsl_vector_get(Psi_T_realImag_func, 0);
	y3Imag = gsl_vector_get(Psi_T_realImag_func, 1);

	FlmkReal = -(0.5)*((y3Real*g_func*csiRel_func/(epsilon_T_noDim_func*Clmk)) + 1.0);

	FlmkImag = -(0.5)*((y3Imag*g_func*csiRel_func/(epsilon_T_noDim_func*Clmk)));

	gsl_vector_set(Flmk_realImag_func, 0, FlmkReal);
	gsl_vector_set(Flmk_realImag_func, 1, FlmkImag);
}//double calculate_Flmk

/*************************************************************************************
 ****** This is taken from Triton, the code that Bart's developed:
 ****** Set up a grid in the true and mean anomaly to calculate the 
 ****** integrands in the definitions of the coefficients clmk, G1lmk,
 ****** G2lmk, G3lmk, and G4lmk  
 **************************************************************************************/
void anomalyGrid(double ecc_func, gsl_vector *trueAn_func, gsl_vector *meanAn_func)
{

	/*************************************************************
	* Set up an array of equidistantly spaced values of the true 
	* anomaly trueAn ranging from 0 to pi
	*************************************************************/
	
	int i = 0;

	double h = 0.0, EE = 0.0, trueAn = 0.0, 
	dtrueAn = PI/static_cast<double>(ngrid-1) ;
	
	
	for (i = 0; i < ngrid; i++) {gsl_vector_set(trueAn_func, i, static_cast<double>(i)*dtrueAn);}

	/*************************************************************
	* Calculate the associated values of the eccentric anomaly EE 
	* and the mean anomaly meanAn. Note that trueAn/2 and EE/2 must    
	* be in the same quadrant!
	*************************************************************/	
	for(i = 0; i < ngrid; i++)
	{
		trueAn = gsl_vector_get(trueAn_func, i);
		
		h = sqrt((1.0-ecc_func)/(1.0+ecc_func));
		EE = 2.0*atan(h*tan(trueAn/2.0));

		if (sin(trueAn)*sin(EE) < 0.0)
			EE = EE - PI;
		
		while(1>0)
		{
			if (EE < 0.0) {EE = EE + 2.0*PI;}
			else if (EE > 2.0*PI) {EE = EE - 2.0*PI;}
			else{break;}
		}
		gsl_vector_set(meanAn_func, i, EE - ecc_func*sin(EE));
	}	

}

/*************************************************************************************
 ****** This routine is used to calculated the Fourier coefficients klmk
 ****** and it is taken from Bart Willems' code Triton
 ****** klmk are calculated using Eq. 53 in Willems et al.(2010) 2010ApJ...713..239W
 **************************************************************************************/
double calculateCoeffs_Klmk(int deg_m_int_func, int deg_k_int_func)
{
	double k_lmk = 0.0;

	if (deg_k_int_func < 0)
		errorMessageAndExit("calculateCoeffs_Klmk", "Unknown size of k");
	else if (deg_k_int_func == 0)
	{
		if (deg_m_int_func < 0){k_lmk = 0.0;}
		else if (deg_m_int_func == 0){k_lmk = 0.5;}			
		else {k_lmk = 1.0;}		
	}
	else
		k_lmk = 1.0;
	
	return k_lmk;

}
/*************************************************************************************
 ****** Calculation of an integral int_a^b f(x)dx using the extended trapezium rule
 **************************************************************************************/
double integrate_with_extended_trapezium_rule (gsl_vector* x, gsl_vector* f, int size) 
{ 
	
	int i = 0;
	double integral = 0.0, xA = 0.0, xB = 0.0, fA = 0.0, fB = 0.0;
	
	integral = 0.0;
	for (i = 1; i < size; i++)
	{
		xA = gsl_vector_get(x, i-1);
		xB = gsl_vector_get(x, i);
		
		fA = gsl_vector_get(f, i-1);
		fB = gsl_vector_get(f, i);
		
		integral = integral + (xB - xA)*(fA + fB)/2.0;
	}
	return integral;
}
/*************************************************************************************
 ****** This routine is used to calculated the tidal timescale da/dt following
 ****** Eq. 54  in Willems et al.(2010) 2010ApJ...713..239W
 **************************************************************************************/

double calculate_dAdT_tides(int deg_l_int_func, int deg_m_int_func, int deg_k_int_func, double ecc_func, double R1_noDim_func, double a_noDim_func, 
							double epsilon_T_noDim_func, double Porb_noDim_func, double M1_noDim_func, double M2_noDim_func, 
							int adiabatic_func, double arg_Flmk_func, double Mod_Flmk_func)
{
		
	double dAdT_tides = 0.0, Glmk_2 = 0.0, Klmk = 0.0, 
	
	deg_l = static_cast<double>(deg_l_int_func);
	
	Klmk = calculateCoeffs_Klmk(deg_m_int_func, deg_k_int_func);
		
	Glmk_2 = calculate_Glmk_2(deg_l_int_func, deg_m_int_func, deg_k_int_func, ecc_func, R1_noDim_func, a_noDim_func);
	
	dAdT_tides = (8.0*PI/Porb_noDim_func)*(M2_noDim_func/M1_noDim_func)*
	a_noDim_func*pow(R1_noDim_func/a_noDim_func, deg_l+3.0)*Klmk*Mod_Flmk_func*sin(arg_Flmk_func)*Glmk_2;
	
	
	return dAdT_tides;
}//double calculate_dAdT

/*************************************************************************************
 ****** This routine is used to calculated the GW driven timescale da/dt following
 ****** Peters (1964)
 **************************************************************************************/

double calculate_dAdT_GR(double a_noDim_func, double M1_noDim_func, double M2_noDim_func, double c_light_noDim_func)
{
	
	double dAdT_GR = 0.0;
	
	dAdT_GR = -(64.0/5.0)*M1_noDim_func * M2_noDim_func*(M1_noDim_func + M2_noDim_func)/(pow(a_noDim_func, 3.0)*pow(c_light_noDim_func, 5.0));
	
	//	cout <<"from Kilic paper = "<< (dAdT_GR*2.0*PI/sqrt(Gsolar*(M1_noDim_func+M2_noDim_func)))*(3.0/2.0)*sqrt(a_noDim_func)<<endl;
	//	cout << "dAdT_GR cazzissimo????? = "<< dAdT_GR << endl;
	
	return dAdT_GR;
}//double calculate_dAdT

 /*************************************************************************************
 ****** This routine is used to calculated the tidal timescale daOmega/dt following
 ****** Eq. 67  in Willems et al.(2010) 2010ApJ...713..239W
 **************************************************************************************/

double calculate_dOmegadT_tides(int deg_l_int_func, int deg_m_int_func, int deg_k_int_func, double ecc_func, double R1_noDim_func, double a_noDim_func, 
								double epsilon_T_noDim_func, double Porb_noDim_func, double M1_noDim_func, double M2_noDim_func, 
								int adiabatic_func, double momInertia_noDim_func, double arg_Flmk_func, double Mod_Flmk_func)
{
	
	const double 
	M1_2 = M1_noDim_func * M1_noDim_func,
	M2_2 = M2_noDim_func * M2_noDim_func,
	Mtot = M1_noDim_func + M2_noDim_func;
	
	double dOmegadT_tides = 0.0, Glmk_4 = 0.0, Klmk = 0.0, Glmk_2 = 0.0, Glmk_3 = 0.0, 
	deg_l = static_cast<double>(deg_l_int_func);
	
	Klmk = calculateCoeffs_Klmk(deg_m_int_func, deg_k_int_func);	
	Glmk_2 = calculate_Glmk_2(deg_l_int_func, deg_m_int_func, deg_k_int_func, ecc_func, R1_noDim_func, a_noDim_func);
	Glmk_3 = calculate_Glmk_3(deg_l_int_func, deg_m_int_func, deg_k_int_func, ecc_func, R1_noDim_func, a_noDim_func);
	
	if(fabs(ecc_func - 0.0) <= 1.0e-10){Glmk_4 = -Glmk_2/2.0;}
	else {Glmk_4 = (ecc_func/sqrt(1.0-ecc_func*ecc_func))*(Glmk_3 - ((1.0-ecc_func*ecc_func)/(2.0*ecc_func))*Glmk_2);}

	dOmegadT_tides = (8.0*PI/Porb_noDim_func)* sqrt(M1_2*M2_2/Mtot)*(M2_noDim_func/M1_noDim_func)*
	(sqrt(a_noDim_func)/momInertia_noDim_func)*pow(R1_noDim_func/a_noDim_func, deg_l + 3.0)*
	Klmk*Mod_Flmk_func*sin(arg_Flmk_func)*Glmk_4;
	
	return dOmegadT_tides;
}//double calculate_dOmegadT_tides

/*************************************************************************************
 ****** This routine is used to calculated the tidal timescale de/dt following
 ****** Eq. 55  in Willems et al.(2010) 2010ApJ...713..239W
 **************************************************************************************/

double calculate_dedT_tides(int deg_l_int_func, int deg_m_int_func, int deg_k_int_func, double ecc_func, double R1_noDim_func, double a_noDim_func, 
							double epsilon_T_noDim_func, double Porb_noDim_func, double M1_noDim_func, double M2_noDim_func, 
							int adiabatic_func, double arg_Flmk_func, double Mod_Flmk_func)

{
	double dedT_tides = 0.0, Glmk_3 = 0.0, Klmk = 0.0, 
	deg_l = static_cast<double>(deg_l_int_func);
	
	Klmk = calculateCoeffs_Klmk(deg_m_int_func, deg_k_int_func);
	
	Glmk_3 = calculate_Glmk_3(deg_l_int_func, deg_m_int_func, deg_k_int_func, ecc_func, R1_noDim_func, a_noDim_func);
	
	dedT_tides = (8.0*PI/Porb_noDim_func)*(M2_noDim_func/M1_noDim_func)*
	pow(R1_noDim_func/a_noDim_func, deg_l+3.0)*Klmk*Mod_Flmk_func*sin(arg_Flmk_func)*Glmk_3;
	
	return dedT_tides;
}//double calculate_dedT_tides


