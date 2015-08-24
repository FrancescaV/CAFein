/***********************************************************************************************************
 !***********************************************************************************************************
 !*****************************									   	  **************************************
 !*****************************		Stellar oscillation code (2011)	  **************************************
 !*****************************									   	  **************************************
 !***********************************************************************************************************
 !*********	
 !*********	    This code has been developed by Francesca Valsecchi in Winter/spring/summer 2011 and 2012
 !*********		to calculate eigenfrequencies and eigenfunctions of WD models and tidally excited stellar
 !*********		pulsations.	
 !*********				
 !*********		The equations for stellar oscillations are the one derived by Wasaburo Unno et al.
 !*********		and explained in the book "Non radial oscillation of stars".
 !*********		The equations are integrated using the Riccati method as described by Takata & Loffler 2004 
 !*********		(2004PASJ...56..645T) with some modifications. Read this paper to understand CAFein.
 !*********		 
 !*********		 Once the tidal potential is included, the equations remain the same, apart for a 
 !*********		 re-definition of the perturbation of the potential and for a different set of 
 !*********		 BC at the surface.
 !*********		 
 !*********		Here is the idea:
 !*********	    
 !*********	    ------------------------- Pure Stellar oscillations --------------------------------
 !*********	    The calculation of eigenfrequencies and eigenfunctions is performed in 2 steps:
 !*********	    
 !*********	    Step 1) Calculation of the eigenfrequencies exactly as described in Takata & Loffler 2004.
 !*********	    The components of the Riccati matrix R are integrated for a set of normalized frequencies w.
 !*********	    To cover the whole star, the integration is performed in two steps: 
 !*********	    - From the surface to a fitting point xFit (Rin)
 !*********	    - From the center to xFit (Rout)
 !*********	    The fitting point is chosen based on the behavior of the brunt vaisala frequency.
 !*********	    The integration is performed with the proper initial conditions (see Wasabuto Unno) at 
 !*********	    the star's boundaries.
 !*********	    During the integration, if the Euclidean norm of the Riccati matrix ||R|| is bigger than 
 !*********	    some user defined limit, then a permutation is applied. 
 !*********	    Specifically, the code calculated the permutation with the minimum ||R||, and re-start 
 !*********	    the integration. 
 !*********	    The normalized frequency w is an eigenfrequency if det (Rout - Rin)_xFit = 0.
 !*********	    
 !*********	    Step 2) Calculation of the eigenfunctions exactly as described in Takata & Loffler 2004.
 !*********	    The calculation of the eigenfunctions require the Riccati components calculated during step 1. 
 !*********	    The integration is performed in two steps: 
 !*********	    - From XFit to the center, 
 !*********	    - From xFit to the surface, 
 !*********	    with proper initial conditions (see Takata & Loffler). The code has the option of either
 !*********	    reading the calculated R and interpolate, or re-integrate R together with the eigenfunctions. 
 !*********     In the second case, some care has to be taken with the integration, because integrating 
 !*********     R from the fitting point to the stars boundaries is numerically unstable. To get around 
 !*********     numerical problems I do the following:
 !*********     During Step 1, every time the code perform a permutation of R I store the information about 
 !*********     that permutation, and the permuted riccati matrix R'. The values that are stored are then used
 !*********     as initial conditions during the calculation of the eigenfunctions. 
 !*********     For this reason, the integration of the eigenfunctions is broken across several sub-intervals 
 !*********     based on how many permutations were performed in step 1.
 !*********	
 !*********	    ------------------------- Tidal Stellar oscillations --------------------------------
 !*********	    The eigenfrequency this time is given in input (fixing spin period and orbital period)
 !*********	    
 !*********	    Step 1) Calculation of the Riccati components as above.
 !*********	    - From the surface to a fitting point xFit (Rin)
 !*********	    - From the center to xFit (Rout)
 !*********	    The fitting point is chosen based on the behavior of the brunt vaisala frequency.
 !*********	    The integration is performed with the proper initial conditions (see Wasabuto Unno and 
 !*********	    Bart’s thesis (Eq. 2.67) for the equation requiring the continuity of the potential) at 
 !*********	    the star's boundaries.
 !*********	    During the integration, if the Euclidean norm of the Riccati matrix ||R|| is bigger than 
 !*********	    some user defined limit, then a permutation is applied (apart very close to the star's surface). 
 !*********	    to reduce noise there.
 !*********	    Specifically, the code calculated the permutation with the minimum ||R||, and re-start 
 !*********	    the integration. 
 !*********	    
 !*********	    Step 2) Calculation of the eigenfunctions exactly as described in Takata & Loffler 2004.
 !*********	    The calculation of the eigenfunctions require the Riccati components calculated during step 1. 
 !*********	    The integration is performed in two steps: 
 !*********	    - From XFit to the center, 
 !*********	    - From xFit to the surface, 
 !*********	    with proper initial conditions (see Takata & Loffler). The riccati components are interpolated
 !*********	
 !*********	    Note, with the inclusion of non adiabatic effects I have to integrate many more equations
 !*********	    accounting for the real and imaginary part. I kept the same structure as above and added more equations
 !*********	    for the non-adiabatic case. 
 !*********	    FEW MORE DETAILS below..
 !*********	    
 !**************************************************************************************************************	    
 !**************************************************************************************************************	    
 !*********	    
 !*********	  1) Original configuration of U and V
 !*********	  
 !*********	  adiabatic, no tides:
 !*********	  
 !*********	  U = y1, y2
 !*********	  V = y3, y4
 !*********	  
 !*********	  adiabatic, tides:
 !*********	  
 !*********	  U = y1, y2, y7
 !*********	  V = y3, y4, y8
 !*********	  
 !*********	  non adiabatic, no tides:
 !*********	  
 !*********	  U = y1r, y2r, y5r, y1i, y2i,y5i
 !*********	  V = y3r, y4r, y6r, y3i, y4i, y6i 
 !*********	  
 !*********	  non adiabatic, tides:
 !*********	  
 !*********	  U = y1r, y2r, y5r, y7r, y1i, y2i, y5i, y7i
 !*********	  V = y3r, y4r, y6r, y8r, y3i, y4i, y6i, y8i 
 !*********	    
 !*********	  Where r and i denote the real and imaginary part.
 !*********	    
 !**************************************************************************************************************	    
 !**************************************************************************************************************	    
 !*********	    
 !*********	  2) DEFAULT INTEGRATORS (integration of V == integration eigenfunctions):
 !*********	    
 !*********	  - Calculation of R, looking for eigenfrequencies: bsimp (in a loop).  
 !*********	    
 !*********	  - Calculation of V:
 !*********	  Adiabatic case:   
 !*********	    integration R = rk4
 !*********	    integration V = rk4
 !*********	    
 !*********	  Non Adiabatic case:   
 !*********	    integration R = bsimp + refinement close to the permutation points
 !*********	    integration V = bsimp
 !*********	    
 !**************************************************************************************************************	    
 !**************************************************************************************************************	    
 !*********	    
 !*********	  3) NOTE ON SIMMETRY OF R FOR NON ADIABATIC CASE.
 !*********	    
 !*********	  By definition U = RV:
 !*********	  
 !*********	  |Ureal|	|R00    R01||Vreal|
 !*********	  |     | = |		   ||	  |
 !*********	  |Uimag|   |R10    R11||Vimag|
 !*********	  
 !*********	  So:
 !*********	  
 !*********	  Ureal = R00 Vreal + R01 Vimag
 !*********	  Uimag = R10 Vreal + R11 Vimag
 !*********	  
 !*********	  With a complex integrator (which I am not using) my Eq U = RV would read as:
 !*********	    
 !*********	  Ureal + i U imag = (Rreal + i Rimag)(Vreal + i Vimag)
 !*********	  
 !*********	  so
 !*********	    
 !*********	  Ureal = Rreal Vreal - Rimag Vimag
 !*********	  Uimag = Rimag Vreal + Rreal Vimag
 !*********	  
 !*********	  so
 !*********	  Rreal = R00 = R11
 !*********	  Rimag = R10 = -R01
 ***********************************************************************************************************
 ***********************************************************************************************************
 ************************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


#include "dma.h"
#include "SteffenInterp.h"
#include "gsl_odeiv_evolve_Riccati.h"
#include "params_integrator_V_nonAd_rkf45.h"
#include "params_integrator_V_nonAd.h"
#include "params_integrator_adiabatic.h"
#include "params_integrator_R_nonAd.h"
#include "matrix_operations.h"
#include "RiccatiMatrices_operations.h"
#include "eigenfunction_operations.h"
#include "IOfiles.h"
#include "initialConditions.h"
#include "readInputParameters.h"
#include "TidalParameters_operations.h"
#include "numericalIntegration.h"
#include "rkf45_state.h"


using namespace std;


int
main (void)
{	
	
	cout << setiosflags(ios_base::scientific) << setprecision(16);
	cerr << setiosflags(ios_base::scientific) << setprecision(16);
	
	//	const double PI = 2.0*acos(0.0), Gsolar =  3.9392318e-07;
	
	const double PI = 2.0*acos(0.0), 
	G =  6.67259e-8, 
	Rsun_toCm  = 6.95660e10, 
	Msun_toG = 1.98855e+33,
	secPerMin = 60.0,
	minsPerYears = 60.0*24.0*365.242199,   //m/h h/d d/y;
	secPerYears = secPerMin*minsPerYears,	   //s/m m/h h/d d/y
	c_cgs = 2.9979245800e10,
	Gsolar = G * Msun_toG/pow(Rsun_toCm, 3.0);
	
	const int Willems2003Test = 0; 
	/***************************************************************	
	 ***** Getting accuracies and various params from an input file.
	 ***** 
	 ***** The parameters to refine the search of an eigenfreq. via
	 ***** the secant method are also given (in case the method is
	 ***** applied). For those parameters detR refers to 
	 ***** det(Rinward - Routward), which has a real and imaginary part 
	 ***** 
	 ***** The parameters to calculate tidally excited oscillations
	 ***** are also given (in case tides are applied). 
	 ***************************************************************/
	int adiabatic = 0, detailedRiccati = 0, secant_omega_ri = 0,
	numStepsOmega_r = 0, interpol_rij = 0, reset_h_eachNsteps = 0,
	keep_hConst_forNsteps = 0, nStart = 0, nEnd = 0, dummy = 0, 
	writeEigenfunctionsSkip = 0, integratorRforV = 0, 
	integratorV = 0, integratorFC = 0, integratorFS = 0,
	nonAdiabatic = 0, numStepsOmega_i = 0, permutationNumber = 0, 
	counterIntState = 0, mesh_computedRij = 0, size_rkf45_state = 0;
	
	double 	deg_l_dbl = 0.0, h_init = 0.0, eps_abs = 0.0, eps_rel = 0.0, LimitToRiccati = 0.0, 
	w_min_r_noDim = 0.0, w_max_r_noDim = 0.0, displaceEfreq = 0.0, eigenv_r = 0.0, eigenv_i = 0.0, 
	h_init_eigenfunc = 0.0, eps_abs_eigenfunc = 0.0, eps_rel_eigenfunc = 0.0, 
	numStepsOmega_r_dbl = 0.0, numStepsOmega_i_dbl = 0.0,
	w_min_i_noDim = 0.0, w_max_i_noDim = 0.0, numStepsPspin_dbl = 0.0, LimitToRiccatiUserDef = 0.0, 
	checkIC_fittingPoint = 0.0, surfaceFor_rkf45 = 0.0, centerFor_rkf45 = 0.0;
	
	/* Secant method parameters*/
	double omegaI_n_1 = 0.0, omegaI_n_2 = 0.0, omegaR_n_2 = 0.0,omegaR_n_1 = 0.0, 
	g_r_n_1 = 0.0, g_r_n_2 = 0.0, g_i_n_1 = 0.0, g_i_n_2 = 0.0, g_r = 0.0, g_i = 0.0, 
	delta_g_r = 0.0, delta_g_i = 0.0, delta_g21_r = 0.0, 
	delta_g21_i = 0.0, delta_w21_r = 0.0, delta_w21_i = 0.0;
	
	/* tides parameters */
	int deg_l_int = 0, deg_m_int = 0, deg_k_int = 0;
	int tidesFlag = 0, numStepsPspin = 0;
	double M2_Msun = 0.0, Porb_minutes = 0.0, Pspin_minutes = 0.0, Pspin_minutes_min = 0.0, 
	Pspin_minutes_max = 0.0, deg_m_dbl = 0, deg_k_dbl = 0, klmk = 0.0, 
	dadt_tides_noDim = 0.0, dedt_tides_noDim = 0.0, dOmegadt_tides_noDim = 0.0, dadt_GR_noDim = 0.0, 
	tidesTimescale_a = 0.0, GRtimescale_a = 0.0, tidesTimescale_e = 0.0, tidesTimescale_Omega = 0.0, 
	Mod_Flmk = 0.0, arg_Flmk = 0.0, ecc = 0.0, Jorb_noDim = 0.0, Jspin_noDim = 0.0, 
	JdotOrb_noDim = 0.0, JdotSpin_noDim = 0.0, 
	Pspin_noDim = 0.0, omegaSpin_noDim = 0.0, omegaOrb_noDim = 0.0, omegaOrb_peri_noDim = 0.0, Porb_noDim = 0.0, 
	M2_noDim = 0.0, Clmk_noDim = 0.0, Glmk_2_noDim = 0.0, Glmk_3_noDim = 0.0;
	
	gsl_complex Flmk, y_complex;
	
	typedef vector<double> Row;
	vector<Row> table_acc_and_dims;
	
	const char *filenameAcc_dims;
	ifstream acc_and_dims;
	filenameAcc_dims = "input/acc_and_dim.dat";		
	
	readInputParams(filenameAcc_dims, adiabatic, detailedRiccati,
					secant_omega_ri, numStepsOmega_r, numStepsOmega_i, interpol_rij, 
					reset_h_eachNsteps, keep_hConst_forNsteps,
					writeEigenfunctionsSkip, integratorRforV, 
					integratorV, integratorFC, integratorFS, 
					nStart, nEnd, deg_l_dbl, h_init, eps_abs, eps_rel, LimitToRiccatiUserDef, w_min_r_noDim, 
					w_max_r_noDim, displaceEfreq, w_min_i_noDim, w_max_i_noDim, eigenv_r, eigenv_i, 
					h_init_eigenfunc, eps_abs_eigenfunc, eps_rel_eigenfunc,
					omegaR_n_2, omegaR_n_1, omegaI_n_2, omegaI_n_1, g_r_n_2, g_r_n_1,
					g_i_n_2, g_i_n_1, g_r, g_i, tidesFlag, numStepsPspin, deg_m_dbl, deg_k_dbl, M2_Msun, 
					Porb_minutes, Pspin_minutes_min, Pspin_minutes_max, deg_l_int, deg_m_int, deg_k_int, ecc);
	
	numStepsOmega_r_dbl	= static_cast<double>(numStepsOmega_r); 
	numStepsOmega_i_dbl	= static_cast<double>(numStepsOmega_i); 
	numStepsPspin_dbl	= static_cast<double>(numStepsPspin); 
	LimitToRiccati = LimitToRiccatiUserDef;
	
	if (adiabatic == 0){nonAdiabatic = 1;}	
	if(secant_omega_ri)
	{		
		cout << "*********************************************" << endl;
		cout << "*********************************************" << endl;
		cout << "Variables needed for the secant method: " << endl;
		cout << "*********************************************" << endl;
		cout << "*********************************************" << endl;
		cout << "root g(w) real =  " << g_r << endl;
		cout << "root g(w) imag =  " << g_i << endl;
		cout << "omega R n-2    =  " << omegaR_n_2 << endl;
		cout << "omega R n-1    =  " << omegaR_n_1 << endl;
		cout << "omega I n-2    =  " << omegaI_n_2 << endl;
		cout << "omega I n-1    =  " << omegaI_n_1 << endl;
		cout << "g(w) real n-2  =  " << g_r_n_2 << endl;
		cout << "g(w) real n-1  =  " << g_r_n_1 << endl;
		cout << "g(w) imag n-2  =  " << g_i_n_2 << endl;
		cout << "g(w) imag n-1  =  " << g_i_n_1 << endl;
		cout << "*********************************************" << endl;
		cout << "*********************************************" << endl;
	}//if(secant_omega_ri)
	if(tidesFlag)
	{
		cout << "*********************************************" << endl;
		cout << "*********************************************" << endl;
		cout << "Tidal potential variables: " << endl;
		cout << "*********************************************" << endl;
		cout << "*********************************************" << endl;
		cout << "Companion mass (Msun) = " << M2_Msun << endl;
		cout << "P orbital (mins) .....= " << Porb_minutes << endl;
		
		if(Willems2003Test)
			cout << "\n Performing a test like Willems+2003 \n" << endl;
		else
		{
			cout << "P spin min (mins) ....= " << Pspin_minutes_min << endl;
			cout << "P spin max (mins) ....= " << Pspin_minutes_max << endl;		
		}
		
		cout << "degree m..............= " << deg_m_dbl << endl;
		cout << "degree k..............= " << deg_k_dbl << endl;
		cout << "eccentricity..........= " << ecc << endl;
		cout << "*********************************************" << endl;
		cout << "*********************************************" << endl;
	}//if(tidesFlag)	
	/*******************************************************
	 *** Some more variables...
	 *******************************************************/
	double omega_r_noDim = 0.0, omega_i_noDim = 0.0, omega_r2_noDim = 0.0, dbl_loopOmega_r = 0.0, 
	dbl_loopOmega_i = 0.0, dbl_loopPspin = 0.0, csiRel = 0.0, csi_final = 0.0, h = 0.0, 
	normRminimum = 0.0, c1_center = 0.0, V_surf = 0.0, delAD_surf = 0.0, 
	det_RoMinusRin = 0.0, centerRel = 0.0, surfaceRel = 0.0, M1_Msun = 0.0, 
	L1_Lsun = 0.0, Teff_K = 0.0, 
	R1_Rsun = 0.0, csiIn_force = 0.0, csiFin_force = 0.0, omega_T_noDim = 0.0, 
	del_surf = 0.0, columnsEigenfuncFile = 0, wBreakUp = 0.0, PbreakUp_min = 0.0,
	yReNormU = 0.0, yImNormU = 0.0, yReNormV = 0.0, yImNormV = 0.0, y8r_atR = 0.0, y8i_atR = 0.0,
	modulus_y8_atR = 0.0, arg_y8_atR = 0.0, modulus_y_atR = 0.0, arg_y_atR = 0.0,
	yReU = 0.0, yImU = 0.0, yReV = 0.0, yImV = 0.0, savedRijIn = 0, savedRijOut = 0, 
	kr2 = 0.0, Cs2 = 0.0, kh2 = 0.0, 
	xi_r_dyn_R = 0.0, xi_r_st_R = 0.0, xi_r_dyn_I = 0.0, xi_r_st_I = 0.0, modXi_r_osc = 0.0, 
	sec_units = 0.0, grams_units = 0.0, cm_units = 0.0, 
	M1_noDim = 0.0, R1_noDim = 0.0, a_noDim = 0.0, epsilon_T_noDim = 0.0, RocheLobe_noDim = 0.0, c_light_noDim = 0.0, 
	rho_surf_noDim = 0.0, g_surf_noDim = 0.0;
	
	long double normR = 0.0, normR_0 = 0.0, normR_1 = 0.0, normR_2 = 0.0;
	
	int permutationApplied = 0, inwardOutward = 0, z = 0, checkIn = 0, checkOut = 0, 
	counterPermutIn = 0, counterPermutOut = 0, counterPermut = 0, mesh_elements_Rin = 0, 
	mesh_elements_Rout = 0, mesh_elements_Rmax = 0, mesh_elements_Rij = 0, mesh_elements_Rpick = 0, startInterpol_rij = 0,endInterpol_rij = 0, 
	sizeRij_calc = 0, 
	PermutationMinNormNumber = 0, counterIntegrationSteps = 0, 
	i = 0, k = 0, j = 0, m = 0, n = 0, loopOmega_r = 0, loopOmega_i = 0, 
	loopPspin = 0, Nint_R = 0, checkMaxNormR = 0, checkMaxNormR_OK = 0, counterIntStepsOmegaR_OmegaI0 = 0, linearityViol = 0, 
	Nint_R_atCheckNorm = 0, Nint_R_tot = 0;
	
	double **Rij_calc_ad, *csi_calc_ad,*mass_model, *csiRel_model, *Vg_model, *Astar_model, 
	*U_model, *c1_model, *rho_model, *g_model, *N2_model, *L2_model, 
	*Vt_model, *V_model, *del_model, *delad_model, *ks_model, *c2_model, *c3_model, 
	*c4_model, *epsAd_model, *dlnLR_dlnr_model, *epsS_model, 
	*P_model, *Gamma1_model, *T_model, *Cp_model, *S_model, *Lr_model, 	
	*deltaRadRegion_model, **fitVgCoeffs, **fitAstarCoeffs, **fitUcoeffs, **fitC1coeffs, 
	**fitVtCoeffs, **fitVCoeffs, **fitDelCoeffs, **fitDelADcoeffs, **fitKsCoeffs, 
	**fitC2coeffs, **fitC3coeffs, **fitC4coeffs, **fitEpsADcoeffs, **fitdlnLR_dlnrCoeffs, 
	**fitEpsScoeffs, *VgFuncs, *AstarFuncs,*UFuncs,*c1Funcs, *VtFuncs, *VFuncs,*delFuncs,
	*delADFuncs, *ksFuncs, *c2Funcs,*c3Funcs, *c4Funcs, *epsADFuncs, *dlnLR_dlnrFuncs, *epsSFuncs,	
	**fitRhoCoeffs, *rhoFuncs, **fit_g_Coeffs, *gFuncs, **fitPcoeffs, *PFuncs, 
	**fitGamma1coeffs, *Gamma1Funcs, **fitN2coeffs, *N2Funcs, 
	**fitTcoeffs, *TFuncs, **fit_Cp_Coeffs, *CpFuncs, 
	**fitScoeffs, *SFuncs, **fit_Lr_Coeffs, *LrFuncs; 
	
	/* For eigenfunctions */
	vector<Row> elements_Rin, elements_Rout, rRel_and_rkf45_stateIn, rRel_and_rkf45_stateOut, dummyTable;
	int flagPrintOut = 0, loopInterpol = 0, IntegrationStepsEigenfunc = 0;		
	
	/*To determine the size of the Riccati space */
	const int MAX_NUMBER_PERMUTATIONS = 1000,
	MAX_NUMBER_INTEGRATIONSTEPS_AD = 500000;	
	int nvar_Riccati = 0, rowsRic = 0, position_Vprime_y_Riccati = 0, 
	nvar_Riccati_V = 0;
	
	rhoFuncs   = dfun(4);		
	gFuncs	   = dfun(4);		
	VgFuncs    = dfun(4);
	AstarFuncs = dfun(4);
	UFuncs     = dfun(4);
	c1Funcs    = dfun(4);
	VtFuncs    = dfun(4);
	VFuncs     = dfun(4);
	delFuncs   = dfun(4);
	delADFuncs = dfun(4);
	ksFuncs    = dfun(4);
	c2Funcs    = dfun(4);
	c3Funcs    = dfun(4);
	c4Funcs    = dfun(4);
	epsADFuncs = dfun(4);
	epsSFuncs  = dfun(4);
	dlnLR_dlnrFuncs = dfun(4);
	PFuncs	   = dfun(4);		
	Gamma1Funcs= dfun(4);		
	N2Funcs	   = dfun(4);		
	TFuncs	   = dfun(4);		
	CpFuncs	   = dfun(4);		
	SFuncs	   = dfun(4);		
	LrFuncs	   = dfun(4);		
	
	for (i=0; i<4; i++)
	{
		rhoFuncs[i]   = 0.0;
		gFuncs[i]     = 0.0;
		VgFuncs[i]    = 0.0;
		AstarFuncs[i] = 0.0;
		UFuncs[i]     = 0.0;
		c1Funcs[i]    = 0.0;		
		VtFuncs[i]    = 0.0;
		VFuncs[i]     = 0.0;
		delFuncs[i]   = 0.0;
		delADFuncs[i] = 0.0;
		ksFuncs[i]    = 0.0;
		c2Funcs[i]    = 0.0;
		c3Funcs[i]    = 0.0;
		c4Funcs[i]    = 0.0;
		epsADFuncs[i] = 0.0;
		epsSFuncs[i]  = 0.0;
		dlnLR_dlnrFuncs[i] = 0.0;
		PFuncs[i]      = 0.0;
		Gamma1Funcs[i] = 0.0;
		N2Funcs[i]	   = 0.0;
		TFuncs[i]      = 0.0;
		CpFuncs[i]     = 0.0;
		SFuncs[i]      = 0.0;
		LrFuncs[i]     = 0.0;
		
	}
	/***********************************************************
	 *****                additional I/O files
	 ***********************************************************/
	const char *filenamePolytrope;
	
	string fileRiccatiOutward, fileRiccatiInward, fileRiccatiAtFit, 
	filePolytrope, fileRiccatiOutward_refined, fileRiccatiInward_refined,
	fileEigenfunctions, fileEigenfunctionsSurface, fileGlobalProperties;
	/* **********************************************************/
	/*				reading address for outputfiles		        */
	/* **********************************************************/
	const char *fileAddressModelName;	
	string addressModel;
	
	ifstream fileAddressModel;
	fileAddressModelName = "input/stellarModelFilesAddress.dat";
	fileAddressModel.open(fileAddressModelName); 	
	
	/*reading line of the input file into a string. The counter
	 is to avoid reading an empty line*/
	
	i=0;
	while(fileAddressModel && i==0)
	{
		getline(fileAddressModel, addressModel);
		i++;
	}	
	fileAddressModel.close();	
	filenamePolytrope = addressModel.c_str();
	
	if(adiabatic)
	{
		if(integratorFC){fileRiccatiOutward	= "output/riccatiMatrices_out.dat";}		
		fileRiccatiInward					= "output/riccatiMatrices_in.dat";
	}
	fileRiccatiAtFit			= "output/riccatiMatrices_csi1.dat";
	fileEigenfunctions			= "output/eigenfunctions.dat";
	fileEigenfunctionsSurface	= "output/eigenfunctionsAtSurface.dat";
	fileGlobalProperties		= "output/globalPropertiesSystem.dat";
	
	cout << "Input model = " << filenamePolytrope << endl;		
	
	ofstream riccatiMcsi1, elmntsRout_f, elmntsRin_f, eigenfunctions, eigenfunctionsSurface, 
	globalPropertiesSystem, rkf45_state_in, rkf45_state_out, dummyFile;
	
	/*open files*/
	riccatiMcsi1.open(fileRiccatiAtFit.c_str());
	eigenfunctions.open(fileEigenfunctions.c_str());
	eigenfunctionsSurface.open(fileEigenfunctionsSurface.c_str());
	globalPropertiesSystem.open(fileGlobalProperties.c_str());
	
	riccatiMcsi1 << setiosflags(ios_base::scientific);
	eigenfunctions << setiosflags(ios_base::scientific);
	eigenfunctionsSurface << setiosflags(ios_base::scientific);
	globalPropertiesSystem << setiosflags(ios_base::scientific);
	/****************************************************************************************
	 ***** 
	 *****		DETERMINING THE SIZE OF THE RICCATI SPACE
	 ***** 
	 ***** - position_Rprime_matPermutations and position_Rprime_y_Riccati are the starting 
	 ***** position for variables in the integrator and in the matrix storing the permutations info.
	 ***** - the number of rows in R is 1/2 the number of integration variables y_i
	 ***** 
	 ****************************************************************************************/
	if(nonAdiabatic)
	{
		if(tidesFlag)
			rowsRic = 8;
		else
			rowsRic = 6;
		
		nvar_Riccati = rowsRic*rowsRic;
		nvar_Riccati_V = rowsRic;
		position_Vprime_y_Riccati = 0;
	}//if(nonAdiabatic)
	else
	{
		if(tidesFlag)
			rowsRic = 3;
		else
			rowsRic = 2;
		
		nvar_Riccati = rowsRic*rowsRic + rowsRic;		
		nvar_Riccati_V = nvar_Riccati;
		position_Vprime_y_Riccati = static_cast<int>(rowsRic*rowsRic);
	
	}//else adiabatic 
	size_t dimensionODE_R = nvar_Riccati;
	size_t dimensionODE_V = nvar_Riccati_V;
	
	if (rowsRic != 2 && rowsRic != 3 && rowsRic != 6 && rowsRic != 8)
		errorMessageAndExit("CAFein.c", "unkown rowsRic size!");
	
	int sizeR = static_cast<int>(rowsRic*rowsRic), 
	rowsT = rowsRic * 2, 
	permutationsInfo = 2 + rowsT + sizeR + rowsT +1,  //"2" accounts for csi_init and csi_fin and "1" b/c of the permutation #
	sizeDetVector = 2,
	position_Rprime_matPermutations = 2 + rowsT, 
	position_Rprime_y_Riccati = 0;
	
	/*Label the output files */
	labelF_Riccati_xFit(riccatiMcsi1, rowsRic);
	
	/*If you want in the file containing the eigenfunctions also the Riccati coefficients, then uncomment
	 the lines right below and comment the one containing "labelF_eigenfunctions_and_unperturbed_model"*/
	
//	labelF_eigenfunctions(eigenfunctions, rowsRic, adiabatic, tidesFlag);
//	labelF_eigenfunctions(eigenfunctionsSurface, rowsRic, adiabatic, tidesFlag);	
	labelF_eigenfunctions_and_unperturbed_model(eigenfunctions, rowsRic, adiabatic, tidesFlag);
	labelF_eigenfunctions_and_unperturbed_model(eigenfunctionsSurface, rowsRic, adiabatic, tidesFlag);	
	
	labelF_globalProperties(globalPropertiesSystem);
	/**********************************************************************
	 ***  
	 ***  BLOCK NEEDED FOR THE CALCULATION OF V:
	 ***  
	 ***  Rij_calc will store the calculated Rij, in case Rij are interpolated
	 ***  during the calculation of V.
	 ***  
	 ***********************************************************************/		
	sizeRij_calc = sizeR;

	size_rkf45_state = static_cast<int>(2*sizeRij_calc+1);

	Rij_calc_ad = dfunc(MAX_NUMBER_INTEGRATIONSTEPS_AD,sizeRij_calc);
	csi_calc_ad = dfun(MAX_NUMBER_INTEGRATIONSTEPS_AD);
	/************************************************************************************
	 ***  
	 ***  RICCATI MATRICES FOR THE EIGENFREQUENCIES (or for calculating Rij)
	 ***  	 
	 ***  Note: matPermInfo_# stores info about the permutations applied during the integration of Rij.
	 ***  The configuration of the unpermuted U and V are explain in Note #1 at the top of this file.	 
	 ***  	 
	 ***  The integration is broken in sub-intervals whenever ||R||>limit.	 
	 ***  Each row of matPermInfo_# contains info about each integration interval. 	 
	 ***  	 
	 ***  The columns are:	 
	 ***  	 
	 ***  csi_i......csi_f....Ttot(i1, i2, i3...)...Rij'(r00, r01....).....Tstep(i1, i2, i3...)....permutation #	 
	 ***  	 
	 ***  Where:
	 ***  -Ttot = permutation of unpermuted variables. 	 
	 ***  -Tstep = permutation between two adjacent integration intervals	 
	 ***  -i1, i2, etc... are the positions of the "1" in T.	 
	 ***  	 
	 ***  The way the indices go is:
	 ***  	 
	 ***  i = 2; i< (2+rowsT)...............................for Ttot
	 ***  i = (2+rowsT); i< (2+rowsT+sizeR).................for Rij' 
	 ***  i = (2+rowsT+sizeR); i< (2+rowsT+sizeR+rowsT).....for Tstep	
	 ***  i = (2+rowsT+sizeR+rowsT) for permutation #
	 ***  	 
	 ***  Note mat_normR: is given in input to the routine that finds the permutation giving the minimum norm. 
	 ***  The elements are:	 
	 ***  	 
	 ***  - the positions where ``1'' is in T
	 ***  - the ||R|| for a given permutation
	 ***  - the permutation number
	 ************************************************************************************/
	gsl_matrix *matT		 =  gsl_matrix_alloc(rowsT, rowsT);
	gsl_matrix *matTapplied  =  gsl_matrix_alloc(rowsT, rowsT);
	gsl_matrix *matM		 =  gsl_matrix_alloc(rowsT, rowsT);
	gsl_matrix *matR    	 =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matdR    	 =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matT00		 =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matT01		 =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matT10		 =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matT11		 =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matRoriginal  =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matA			   =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matB			   =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matC	           =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matCR	           =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matD		       =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matR_in_xFit	   =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matR_out_xFit      =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matR_diff          =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matR_diff_SVD      =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matV_SVD           =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matR_diff_rr       =  gsl_matrix_alloc(static_cast<int>(rowsRic/2), static_cast<int>(rowsRic/2));
	gsl_matrix *matR_diff_ii       =  gsl_matrix_alloc(static_cast<int>(rowsRic/2), static_cast<int>(rowsRic/2));
	gsl_matrix *matR_diff_ri       =  gsl_matrix_alloc(static_cast<int>(rowsRic/2), static_cast<int>(rowsRic/2));
	gsl_matrix *matR_diff_ir       =  gsl_matrix_alloc(static_cast<int>(rowsRic/2), static_cast<int>(rowsRic/2));
	gsl_matrix *matDummySizeT	    =  gsl_matrix_alloc(rowsT, rowsT);
	gsl_matrix *matDummySizeT_2     =  gsl_matrix_alloc(rowsT, rowsT);
	gsl_matrix *matDummySizeR	    =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matDummySizeR_2     =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matDummySizeR_3     =  gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matRelmnts_in       =  gsl_matrix_alloc(3, sizeR + 2);
	gsl_permutation *permDummySizeR = gsl_permutation_alloc(rowsRic);
	gsl_permutation *permDummySizeT = gsl_permutation_alloc(rowsT);
	gsl_matrix *matPermInfo_in		= gsl_matrix_alloc(MAX_NUMBER_PERMUTATIONS, permutationsInfo);
	gsl_matrix *matPermInfo_out		= gsl_matrix_alloc(MAX_NUMBER_PERMUTATIONS, permutationsInfo);
	gsl_matrix *matPermInfo			= gsl_matrix_alloc(MAX_NUMBER_PERMUTATIONS, permutationsInfo);	
	gsl_vector *vecPermutIndices	= gsl_vector_alloc(rowsT);
	gsl_vector *vec_y_Riccati_all	= gsl_vector_alloc(sizeR);	
	gsl_vector *vec_y_Riccati_all_test	= gsl_vector_alloc(sizeR);	
	gsl_vector *vecPermutInit_BC	= gsl_vector_alloc(rowsT);
 	gsl_vector *vecDummySizeT		= gsl_vector_alloc(rowsT);
	gsl_vector *vecDummySizeT2		= gsl_vector_alloc(rowsT*rowsT);
	gsl_vector *vecDummySizeR		= gsl_vector_alloc(rowsRic);
	gsl_vector *vecDummySizeR_3		= gsl_vector_alloc(rowsRic);
	gsl_vector *vecDummySizeR2		= gsl_vector_alloc(sizeR);	
	gsl_vector *vecDetDiff_r_i		= gsl_vector_alloc(sizeDetVector);	
	gsl_vector *vecS_SVD	     	= gsl_vector_alloc(rowsRic);	
	gsl_vector *work_SVD	     	= gsl_vector_alloc(rowsRic);	
	/********************************************************
	 *****  Riccati matrices for the eigenfunctions...
	 *****        
	 ***** Note: matRprime, vecVprime, and vecUprime 
	 ***** are used to make the code and procedure more readable.
	 *********************************************************/
	gsl_vector *vecU         = gsl_vector_alloc(rowsRic);
	gsl_vector *vecV         = gsl_vector_alloc(rowsRic);
	gsl_vector *vecUoriginal = gsl_vector_alloc(rowsRic);
	gsl_vector *vecVoriginal = gsl_vector_alloc(rowsRic);
	gsl_matrix *matRinit	 = gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_matrix *matRprime    = gsl_matrix_alloc(rowsRic, rowsRic);
	gsl_vector *vecUprime	 = gsl_vector_alloc(rowsRic);
	gsl_vector *vecVprime    = gsl_vector_alloc(rowsRic);
	gsl_matrix *matTstep     = gsl_matrix_alloc(rowsT, rowsT);	
	gsl_matrix *matTtot      = gsl_matrix_alloc(rowsT, rowsT);
	gsl_vector *vecDummySizeR_1 = gsl_vector_alloc(rowsRic);
	gsl_vector *vecDummySizeR_2 = gsl_vector_alloc(rowsRic);	
	gsl_vector *vecUoriginal_normalized = gsl_vector_alloc(rowsRic);
	gsl_vector *vecVoriginal_normalized = gsl_vector_alloc(rowsRic);
	gsl_vector *Flmk_realImag  = gsl_vector_alloc(2);	
	gsl_vector *Psi_T_realImag = gsl_vector_alloc(2);	

	gsl_vector *vec_rkf45_state_in = NULL;
	gsl_matrix *mat_rkf45_state_in = NULL;	
	gsl_vector *vec_rkf45_state_out = NULL;

	if(integratorRforV == 3)
	{
		vec_rkf45_state_in  = gsl_vector_alloc(size_rkf45_state);
		mat_rkf45_state_in	= gsl_matrix_alloc(3, size_rkf45_state);	

		if(integratorFC)
			vec_rkf45_state_out	= gsl_vector_alloc(size_rkf45_state);
	}//if(integratorRforV == 3)
	/*******************************************************
	 *** Riccati integrator stuff...
	 *******************************************************/
	const gsl_odeiv_step_type * T_Riccati; 
	const gsl_odeiv_step_type * T_eigenfunc; 
	
	if(integratorRforV == 1){T_Riccati = gsl_odeiv_step_rk4;}
	else if(integratorRforV == 2){T_Riccati = gsl_odeiv_step_bsimp;}
	else {T_Riccati = gsl_odeiv_step_rkf45;}
	
	if(integratorV == 1){T_eigenfunc = gsl_odeiv_step_rk4;}
	else if (integratorV == 2){T_eigenfunc = gsl_odeiv_step_bsimp;}
	else{T_eigenfunc = gsl_odeiv_step_rkf45;}
	
	gsl_odeiv_step * s_Riccati    = gsl_odeiv_step_alloc (T_Riccati, dimensionODE_R);
	gsl_odeiv_control * c_Riccati = gsl_odeiv_control_y_new (eps_abs,eps_rel);
	gsl_odeiv_evolve * e_Riccati  = gsl_odeiv_evolve_alloc (dimensionODE_R);
	
	gsl_odeiv_step * s_eigenfunc    = gsl_odeiv_step_alloc (T_eigenfunc, dimensionODE_V);
	gsl_odeiv_control * c_eigenfunc = gsl_odeiv_control_y_new (eps_abs_eigenfunc,eps_rel_eigenfunc);
	gsl_odeiv_evolve * e_eigenfunc  = gsl_odeiv_evolve_alloc (dimensionODE_V);
	
	cout << "*************************************" << endl;
	cout << "Stepping method for Riccati = " << gsl_odeiv_step_name (s_Riccati)	<< endl;
	cout << "Stepping method for eigenfunctions = " << gsl_odeiv_step_name (s_eigenfunc)	<< endl;
	cout << "*************************************" << endl;
	
	double y_Riccati[nvar_Riccati], y_Riccati_V[nvar_Riccati_V];
	for (i = 0; i < nvar_Riccati; i++){y_Riccati[i] = 0.0;}
	for (i = 0; i < nvar_Riccati_V; i++){y_Riccati_V[i] = 0.0;}
	
	double *y_Riccati_all;
	y_Riccati_all = dfun(sizeR);
	int status;
	/************************************************************************************
	 ***  
	 ***  READING THE STELLAR MODEL:
	 ***  	 
	 *** The weird looping indices account for the fact that the input file has labels.
	 ***  	 
	 **********************************************************************************/	
	vector<Row> tablePoly;
	int meshModel = readInputFile(filenamePolytrope, tablePoly);
	meshModel--;
	
	csiRel_model = dfun(meshModel);
	mass_model   = dfun(meshModel);
	Vg_model 	 = dfun(meshModel);
	Astar_model  = dfun(meshModel);
	U_model      = dfun(meshModel);
	c1_model     = dfun(meshModel);
	rho_model    = dfun(meshModel);
	g_model      = dfun(meshModel);
	N2_model	 = dfun(meshModel);
	L2_model	 = dfun(meshModel);
	Vt_model	 = dfun(meshModel);
	V_model		 = dfun(meshModel);
	del_model	 = dfun(meshModel);
	delad_model	 = dfun(meshModel);
	ks_model	 = dfun(meshModel);
	c2_model	 = dfun(meshModel);
	c3_model	 = dfun(meshModel);
	c4_model	 = dfun(meshModel);
	epsAd_model	 = dfun(meshModel);
	dlnLR_dlnr_model    = dfun(meshModel);
	epsS_model			= dfun(meshModel); 
	deltaRadRegion_model= dfun(meshModel); 
	P_model             = dfun(meshModel);
	Gamma1_model        = dfun(meshModel);
	T_model             = dfun(meshModel);
	Cp_model            = dfun(meshModel);
	S_model             = dfun(meshModel);
	Lr_model            = dfun(meshModel);
	
	fitRhoCoeffs      = dfunc(4,meshModel);
	fit_g_Coeffs      = dfunc(4,meshModel);
	fitVgCoeffs       = dfunc(4,meshModel);
	fitAstarCoeffs    = dfunc(4,meshModel);
	fitUcoeffs        = dfunc(4,meshModel);
	fitC1coeffs       = dfunc(4,meshModel);
	fitVtCoeffs       = dfunc(4,meshModel);
	fitVCoeffs		  = dfunc(4,meshModel);
	fitDelCoeffs      = dfunc(4,meshModel);
	fitDelADcoeffs    = dfunc(4,meshModel);
	fitKsCoeffs       = dfunc(4,meshModel);
	fitC2coeffs       = dfunc(4,meshModel);
	fitC3coeffs       = dfunc(4,meshModel);
	fitC4coeffs       = dfunc(4,meshModel);
	fitEpsADcoeffs    = dfunc(4,meshModel);
	fitEpsScoeffs     = dfunc(4,meshModel);
	fitdlnLR_dlnrCoeffs = dfunc(4,meshModel);
	fitPcoeffs        = dfunc(4,meshModel);
	fitGamma1coeffs   = dfunc(4,meshModel);
	fitN2coeffs       = dfunc(4,meshModel);
	fitTcoeffs		  = dfunc(4,meshModel);
	fit_Cp_Coeffs	  = dfunc(4,meshModel);
	fitScoeffs		  = dfunc(4,meshModel);
	fit_Lr_Coeffs	  = dfunc(4,meshModel);
	
	
	M1_Msun = tablePoly[1][0];
	R1_Rsun = tablePoly[1][1];
	L1_Lsun = tablePoly[1][2];
	Teff_K  = tablePoly[1][6];
	
	sec_units   = tablePoly[1][7];
	grams_units = tablePoly[1][8];
	cm_units    = tablePoly[1][9];
	
	M1_noDim = M1_Msun * Msun_toG/grams_units;
	R1_noDim = R1_Rsun * Rsun_toCm/cm_units;
	
	/*breakup velocity, from equating the force of gravity at the surface F = Mg = GM^2/R^2 to the centrifugal force f = ma = mv^2/R = mw^2R^2/R = mw^2R*/
	wBreakUp = sqrt(Gsolar*M1_Msun/pow(R1_Rsun, 3.0));
	PbreakUp_min = 2.0*PI/(wBreakUp*60.0);
	
	j = 0;
	for (i=4; i<=meshModel+3; i++)
	{
		mass_model[j]     = tablePoly[i][0];
		rho_model[j]	  = pow(10.0,tablePoly[i][1]);
		P_model[j]        = pow(10.0,tablePoly[i][2]);
		T_model[j]		  = pow(10.0,tablePoly[i][3]);
		csiRel_model[j]   = tablePoly[i][4];
		Vg_model[j]		  = tablePoly[i][5];
		Astar_model[j] 	  = tablePoly[i][6];
		U_model[j]        = tablePoly[i][7];
		c1_model[j]       = tablePoly[i][8];
		N2_model[j]		  = tablePoly[i][9];
		L2_model[j]		  = tablePoly[i][10];
		g_model[j]		  = tablePoly[i][11];
		Lr_model[j]       = tablePoly[i][12];
		Cp_model[j]       = tablePoly[i][16];
		deltaRadRegion_model[j]= tablePoly[i][27] - tablePoly[i][26];		
		Gamma1_model[j]   = tablePoly[i][28];
		S_model[j]        = tablePoly[i][29];
		
		if(nonAdiabatic)
		{
			V_model[j]			= tablePoly[i][13];
			delad_model[j]		= tablePoly[i][14];
			del_model[j]		= tablePoly[i][15];
			Vt_model[j]			= tablePoly[i][17];
			ks_model[j]			= tablePoly[i][18];
			epsAd_model[j]		= tablePoly[i][19];
			epsS_model[j]		= tablePoly[i][20];
			c2_model[j]			= tablePoly[i][21];
			c3_model[j]			= tablePoly[i][22];
			c4_model[j]			= tablePoly[i][23];
			dlnLR_dlnr_model[j] = tablePoly[i][24];
		}//if(nonAdiabatic)
		j++;
	}
	tablePoly.clear();
	
	rho_surf_noDim = rho_model[meshModel-1];
	g_surf_noDim = g_model[meshModel-1];
	
	cout << "g_surf_noDim ........= " << g_surf_noDim << endl;	
	cout << "log(rho_surf_noDim)..= " << log10(rho_surf_noDim) << endl;	
	
	/************************************************************************************
	 ***  Calculating the moment of intertia (int_0^m r^2dm) using Gill & Miller.
	 ***  (needed to calculated the spin tidal timescale)
	 ***  	 
	 *** With my dimensionless quantities: I = 4 PI int_0^1 rho r^4 dr
	 *** I thus determined is dimensionless
	 **********************************************************************************/	
	double *inertia_x, *inertia_f, momInertia_noDim = 0.0, err_momInertia = 0.0;
	
	inertia_x = dfun(meshModel);
	inertia_f = dfun(meshModel);
	
	for(i = 0; i<meshModel; i++)
	{
		inertia_x[i] = csiRel_model[i];
		
		inertia_f[i] = 4.0*PI*pow(csiRel_model[i], 4.0)*rho_model[i];
	}	
	
	momInertia_noDim = FourPt_DoubleStar(inertia_x, inertia_f, meshModel, 0, meshModel-1, err_momInertia);
	cout << "I dimensionless..... = " << momInertia_noDim << "-- Integral error = " << err_momInertia << endl;
	
	free(inertia_x);
	free(inertia_f);
	/************************************************************************************/
	/************************************************************************************/	
	
	/*to avoid problems with Steffen I move the star's boundaries by a bit*/
	centerRel = csiRel_model[0]+(csiRel_model[1] - csiRel_model[0])/displaceEfreq;
	
	surfaceRel = csiRel_model[meshModel-1]-(csiRel_model[meshModel-1] - csiRel_model[meshModel-2])/displaceEfreq;

	cout << "I(Msun Rsun^2) ......= " << momInertia_noDim*grams_units*pow(cm_units,2.0)/(Msun_toG*pow(Rsun_toCm, 2.0)) << endl;
	cout << "I(cgs) ..............= " << momInertia_noDim*grams_units*pow(cm_units,2.0) << endl;
	cout << "centerRel ...........= " << centerRel << endl;	
	cout << "surfaceRel ..........= " << surfaceRel << endl;	
	cout << "Total mass (Msun) ...= " << M1_Msun << endl;	
	cout << "Total mass (noDim) ..= " << M1_noDim << endl;	
	cout << "Total radius (Rsun) .= " << R1_Rsun << endl;	
	cout << "Total radius (noDim).= " << R1_noDim << endl;	
	cout << "Break up spin (min) .= " << PbreakUp_min << endl;	
	cout << "degree deg_l_dbl.....= " << deg_l_dbl << endl;
	cout << "rowsT ...............= " << rowsT << endl;	
	cout << "rowsR ...............= " << rowsRic << endl;	
	cout << "sizeR ...............= " << sizeR << endl;	
	cout << "detailed Riccati.....= " << detailedRiccati << endl;
	cout << "Limit to Riccati.....= " << LimitToRiccati << endl;	
	/************************************************************************************
	 ***  Calculating fitting coeffs for model's variables needed by riccati integrator.
	 **********************************************************************************/	
	calcSteffenInterp(meshModel, csiRel_model, Vg_model, fitVgCoeffs);
	calcSteffenInterp(meshModel, csiRel_model, Astar_model, fitAstarCoeffs);
	calcSteffenInterp(meshModel, csiRel_model, U_model, fitUcoeffs);
	calcSteffenInterp(meshModel, csiRel_model, c1_model, fitC1coeffs);
	calcSteffenInterp(meshModel, csiRel_model, rho_model, fitRhoCoeffs);
	calcSteffenInterp(meshModel, csiRel_model, g_model, fit_g_Coeffs);
	calcSteffenInterp(meshModel, csiRel_model, P_model, fitPcoeffs);
	calcSteffenInterp(meshModel, csiRel_model, Gamma1_model, fitGamma1coeffs);
	calcSteffenInterp(meshModel, csiRel_model, N2_model, fitN2coeffs);
	calcSteffenInterp(meshModel, csiRel_model, T_model, fitTcoeffs);
	calcSteffenInterp(meshModel, csiRel_model, Cp_model, fit_Cp_Coeffs);
	calcSteffenInterp(meshModel, csiRel_model, S_model, fitScoeffs);
	calcSteffenInterp(meshModel, csiRel_model, Lr_model, fit_Lr_Coeffs);
	
	if(nonAdiabatic)
	{
		calcSteffenInterp(meshModel, csiRel_model, V_model, fitVCoeffs);
		calcSteffenInterp(meshModel, csiRel_model, delad_model, fitDelADcoeffs);
		calcSteffenInterp(meshModel, csiRel_model, del_model, fitDelCoeffs);
		calcSteffenInterp(meshModel, csiRel_model, Vt_model, fitVtCoeffs);
		calcSteffenInterp(meshModel, csiRel_model, ks_model, fitKsCoeffs);
		calcSteffenInterp(meshModel, csiRel_model, epsAd_model, fitEpsADcoeffs);
		calcSteffenInterp(meshModel, csiRel_model, epsS_model, fitEpsScoeffs);
		calcSteffenInterp(meshModel, csiRel_model, c2_model, fitC2coeffs);
		calcSteffenInterp(meshModel, csiRel_model, c3_model, fitC3coeffs);
		calcSteffenInterp(meshModel, csiRel_model, c4_model, fitC4coeffs);
		calcSteffenInterp(meshModel, csiRel_model, dlnLR_dlnr_model, fitdlnLR_dlnrCoeffs);
		evalSteffenInterp(meshModel,csiRel_model,fitVCoeffs,surfaceRel,VFuncs);
		evalSteffenInterp(meshModel,csiRel_model,fitDelADcoeffs,surfaceRel,delADFuncs);
		evalSteffenInterp(meshModel,csiRel_model,fitDelCoeffs,surfaceRel,delFuncs);
	}//if(nonAdiabatic)
	
	evalSteffenInterp(meshModel,csiRel_model,fitC1coeffs,centerRel,c1Funcs);
	
	c1_center   = c1Funcs[0];
	V_surf		= VFuncs[0];
	delAD_surf  = delADFuncs[0];
	del_surf	= delFuncs[0];
	
	evalSteffenInterp(meshModel,csiRel_model,fit_g_Coeffs,surfaceRel,gFuncs);
	evalSteffenInterp(meshModel,csiRel_model,fitRhoCoeffs,surfaceRel,rhoFuncs);

	g_surf_noDim = gFuncs[0];
	rho_surf_noDim = rhoFuncs[0];
	/************************************************************************************
	 ***  Structure containing the parameters for the integrator:
	 ***  OnwardOutward is any # >1 and is used only during the integration of the eigenfunctions.
	 **********************************************************************************/	
	inwardOutward = 2;
	
	params_integrator_adiabatic_struct params_eigenfreq_adiabatic = {csiIn_force, csiFin_force, 
		rowsRic, meshModel, omega_r_noDim, omega_i_noDim, deg_l_dbl, VgFuncs, AstarFuncs, UFuncs, c1Funcs, 
		VFuncs, delADFuncs, delFuncs, VtFuncs, ksFuncs, epsADFuncs, epsSFuncs, c2Funcs, 
		c3Funcs, c4Funcs, dlnLR_dlnrFuncs, fitVgCoeffs, fitAstarCoeffs, fitUcoeffs, 
		fitC1coeffs, csiRel_model, matA, matB, matC, matD, matCR, matR, matdR, matT, matM, matDummySizeR, 
		matDummySizeR_2, matDummySizeT, matDummySizeT_2, vecDummySizeR2, vecDummySizeT2, vecDummySizeR, 
		vecDummySizeR_3, permDummySizeT, inwardOutward, interpol_rij, Rij_calc_ad, csi_calc_ad, 
		startInterpol_rij, endInterpol_rij, nonAdiabatic, tidesFlag};
	
	params_integrator_adiabatic_struct params_eigenfunc_adiabatic = {csiIn_force, csiFin_force, rowsRic, meshModel, omega_r_noDim, 
		omega_i_noDim, deg_l_dbl, VgFuncs, AstarFuncs, UFuncs, c1Funcs, VFuncs, delADFuncs, delFuncs, VtFuncs, ksFuncs, 
		epsADFuncs, epsSFuncs, c2Funcs, c3Funcs, c4Funcs, dlnLR_dlnrFuncs, fitVgCoeffs, fitAstarCoeffs, fitUcoeffs, 
		fitC1coeffs, csiRel_model, matA, matB, 
		matC, matD, matCR, matR, matdR, matTtot, matM, matDummySizeR, matDummySizeR_2, matDummySizeT, 
		matDummySizeT_2, vecDummySizeR2, vecDummySizeT2, vecDummySizeR, vecDummySizeR_3, permDummySizeT,
		inwardOutward, interpol_rij, Rij_calc_ad, csi_calc_ad, startInterpol_rij, endInterpol_rij, nonAdiabatic, tidesFlag};
	
	params_integrator_R_nonAd_struct params_eigenfreq_nonAdiabatic = {rowsRic, meshModel, omega_r_noDim, omega_i_noDim, deg_l_dbl, VgFuncs, AstarFuncs, UFuncs, c1Funcs, 
		VFuncs, delADFuncs, delFuncs, VtFuncs, ksFuncs, epsADFuncs, epsSFuncs, c2Funcs, 
		c3Funcs, c4Funcs, dlnLR_dlnrFuncs, fitVgCoeffs, fitAstarCoeffs, fitUcoeffs, 
		fitC1coeffs, fitVCoeffs, fitDelADcoeffs, fitDelCoeffs, fitVtCoeffs, fitKsCoeffs, 
		fitEpsADcoeffs, fitEpsScoeffs, fitC2coeffs, fitC3coeffs, fitC4coeffs, fitdlnLR_dlnrCoeffs, 
		csiRel_model, matA, matB, matC, matD, matCR,matR, matdR, matT, matM, matDummySizeR, 
		matDummySizeR_2, matDummySizeT, matDummySizeT_2, vecDummySizeR2, vecDummySizeT2, vecDummySizeR, 
		vecDummySizeR_3, permDummySizeT, inwardOutward, nonAdiabatic, tidesFlag};
	
	/************************************************************************************
	 ***  If adiabatic I use the integrator that handles both R and V at the same time.
	 ***  In non adiabatic conditions I use two separate integrators
	 **********************************************************************************/	
	gsl_odeiv_system sys_Riccati_R_adiabatic		   = {func_Riccati_R, jac_Riccati_R, dimensionODE_R, &params_eigenfreq_adiabatic};  	
	gsl_odeiv_system sys_Riccati_R_eigenfunc_adiabatic = {func_Riccati_R, jac_Riccati_R, dimensionODE_V, &params_eigenfunc_adiabatic};  	
	
	gsl_odeiv_system sys_Riccati_R_nonAdiabatic		   = {func_Riccati_Ronly, jac_Riccati_Ronly, dimensionODE_R, &params_eigenfreq_nonAdiabatic};
	gsl_odeiv_system sys_Riccati_R_toUse;
	/************************************************************************************
	 ************************************************************************************
	 ************    
	 ************  Calculating the elements of R as described in Takata Loffler 2004 
	 ************  (2004PASJ...56..645T)
	 ************           
	 ************  The configuration of U and V is explained in Note #1 at the top of this file.        
	 ************  	 
	 ************  If non adiabatic, the integrator variables y[i] contain Rij.
	 ************  	 
	 ************  If adiabatic, the integrator variables y[i] contain Rij and V, for instance:
	 ************  	y[0] = r00, y[1] = r01, y[2] = r02
	 ************  	y[3] = r10, y[4] = r11, y[5] = r12, 
	 ************  	y[6] = r20, y[7] = r21, y[8] = r22,
	 ************  	y[9] = v0, y[10] = v1, y[11] = v2
	 ************  	
	 ************  	During the first step I only calculate Rij ==> v_i = 0 as Initial conditions
	 ************  	
	 ************************************************************************************
	 **********************************************************************************/	
	for(loopPspin = 0; loopPspin < numStepsPspin; loopPspin++)
	{
		inwardOutward = 2; 
		/*******************************************************
		 *** opening the file that will store Rij 
		 *** (needed if interpolating R during integration of V)
		 *******************************************************/			
		elmntsRin_f.open(fileRiccatiInward.c_str());
		elmntsRin_f  << setiosflags(ios_base::scientific);
		labelF_Riccati_inOut(elmntsRin_f, rowsRic, tidesFlag);

		if(integratorFC)
		{
			elmntsRout_f.open(fileRiccatiOutward.c_str());
			elmntsRout_f << setiosflags(ios_base::scientific);	
			labelF_Riccati_inOut(elmntsRout_f, rowsRic, tidesFlag);		
		}				
		/*******************************************************
		 *** Calculating the tidal frequency if tides are applied
		 *** 
		 *** sigma = kn + m Omega (calculated below in seconds)
		 *** 
		 *** n is the mean motions (2pi/Porb)
		 *** Omega is the spin frequency (2pi/Pspin)
		 *******************************************************/			
		dbl_loopPspin = static_cast<double>(loopPspin);
		
		Pspin_minutes = (dbl_loopPspin * Pspin_minutes_max + (numStepsPspin_dbl - dbl_loopPspin) * Pspin_minutes_min)/numStepsPspin_dbl;
		
		Pspin_noDim     = Pspin_minutes * secPerMin/sec_units;
		Porb_noDim      = Porb_minutes * secPerMin/sec_units;
		omegaSpin_noDim = 2.0*PI/Pspin_noDim;
		omegaOrb_noDim  = 2.0*PI/Porb_noDim;
		
		cout << "Pspin_noDim ...........= " << Pspin_noDim << endl;
		cout << "Pspin (min) ...........= " << Pspin_noDim * sec_units/secPerMin << endl;
		cout << "Porb_noDim ............= " << Porb_noDim << endl;
		cout << "Porb (min) ............= " << Porb_noDim * sec_units/secPerMin << endl;
		cout << "omegaSpin_noDim .......= " << omegaSpin_noDim << endl;
		cout << "omegaOrb_noDim ........= " << omegaOrb_noDim << endl;
		/************************************************************************************
		 ***  Calculating of the orbital angular velocity at periastron
		 **********************************************************************************/			
		omegaOrb_peri_noDim = (sqrt(1.0 + ecc)/pow(1.0 - ecc, 3.0/2.0))*(2.0*PI/Porb_noDim);		
		
		if(Willems2003Test)
		{
			omegaSpin_noDim = omegaOrb_peri_noDim/2.0;
			Pspin_noDim     = 2.0 * PI/omegaSpin_noDim;
		}
		
		if(tidesFlag)
		{
			omega_T_noDim = deg_k_dbl * omegaOrb_noDim + deg_m_dbl * omegaSpin_noDim;
			
			M2_noDim = M2_Msun * Msun_toG/grams_units;
			
			a_noDim  = pow(Porb_noDim/(2.0*PI), 2.0/3.0)*pow(M1_noDim + M2_noDim, 1.0/3.0);
			
			epsilon_T_noDim = (M2_noDim/M1_noDim)*pow(R1_noDim/a_noDim, 3.0);			
			
			RocheLobe_noDim = a_noDim*(1.0 - ecc)*0.49*pow(M1_noDim/M2_noDim, 2.0/3.0)/(0.6*pow(M1_noDim/M2_noDim, 2.0/3.0)+log(1.0+pow(M1_noDim/M2_noDim, 1.0/3.0)));
			
			klmk = calculateCoeffs_Klmk(deg_m_int, deg_k_int);
			Clmk_noDim = calculateFourierCoeffs_Clmk(deg_l_int, deg_m_int, deg_k_int, ecc, R1_noDim, a_noDim);
			Glmk_2_noDim = calculate_Glmk_2(deg_l_int, deg_m_int, deg_k_int, ecc, R1_noDim, a_noDim);
			Glmk_3_noDim = calculate_Glmk_3(deg_l_int, deg_m_int, deg_k_int, ecc, R1_noDim, a_noDim);
			
			cout << "omegaT (noDim).........= " << omega_T_noDim << endl;
			cout << "omegaOrb_peri (1/sec)..= " << omegaOrb_peri_noDim/sec_units << endl;
			cout << "omegaSpin (1/sec)......= " << omegaSpin_noDim/sec_units << endl;
			cout << "a (Rsun) ..............= " << a_noDim * cm_units/Rsun_toCm << endl;
			cout << "rocheLobe_peri (Rsun)..= " << RocheLobe_noDim * cm_units/Rsun_toCm << endl;
			cout << "epsilon_T (noDim) .....= " << epsilon_T_noDim << endl;
			cout << "Clkm ..................= " << Clmk_noDim << endl;
			cout << "Glkm_2.................= " << Glmk_2_noDim << endl;
			cout << "Glkm_3.................= " << Glmk_3_noDim << endl;			
			cout << "klmk...................= " << klmk << endl;
			/* and force the loop on omega_r_noDim at this value*/
			w_min_r_noDim = omega_T_noDim;
			w_max_r_noDim = omega_T_noDim;
		}//if(tidesFlag)
		
		/*******************************************************************************
		 ***  Looping on the REAL Omega
		 *******************************************************************************/			
		counterIntStepsOmegaR_OmegaI0 = 0;
		for (loopOmega_r = 0; loopOmega_r < numStepsOmega_r; loopOmega_r++)
		{
			permutationNumber = 0;
			dbl_loopOmega_r = static_cast<double>(loopOmega_r);
			omega_r_noDim = (dbl_loopOmega_r*w_max_r_noDim + (numStepsOmega_r_dbl-dbl_loopOmega_r)*w_min_r_noDim)/numStepsOmega_r_dbl;		
			omega_r2_noDim = omega_r_noDim*omega_r_noDim;
			
			/***************************************************************************
			 *** if the SECANT METHOD for finding a root is applied, 
			 *** then calculating omega_r_noDim and omega_i_noDim from the secant input parameters
			 ***************************************************************************/			
			if(secant_omega_ri)
			{
				numStepsOmega_i = 1;
				numStepsOmega_r = 1;
				
				delta_g_r = g_r - g_r_n_2;
				delta_g_i = g_i - g_i_n_2;
				delta_g21_r = g_r_n_2 - g_r_n_1;
				delta_g21_i = g_i_n_2 - g_i_n_1;
				delta_w21_r = omegaR_n_2 - omegaR_n_1;
				delta_w21_i = omegaI_n_2 - omegaI_n_1;			
				
				/*omega real*/
				omega_r_noDim  = omegaR_n_2 + (1.0/(delta_g21_r*delta_g21_r+delta_g21_i*delta_g21_i))*
				((delta_g_r*delta_w21_r - delta_g_i*delta_w21_i)*delta_g21_r + 
				 (delta_g_r*delta_w21_i + delta_g_i*delta_w21_r)*delta_g21_i);
				
				/*omega imaginary*/
				omega_i_noDim  = omegaI_n_2 + (1.0/(delta_g21_r*delta_g21_r+delta_g21_i*delta_g21_i))*
				((delta_g_r*delta_w21_i + delta_g_i*delta_w21_r)*delta_g21_r - 
				 (delta_g_r*delta_w21_r - delta_g_i*delta_w21_i)*delta_g21_i);
				
				/*omega real squared used to find the fitting point*/
				omega_r2_noDim = omega_r_noDim*omega_r_noDim;
				
				if (fabs(omega_i_noDim) < 1.0e-15){omega_i_noDim = 0;}
				
				/*to force the loop in omega_i_noDim between the values just calculated*/
				w_max_i_noDim = omega_i_noDim;
				w_min_i_noDim = omega_i_noDim;
				
				/*re-setting the initial time step for the integrator*/
				if(n==0){h = -h_init;}
				else{h = h_init;}
			}
			/*****************************************************************************************
			 ***  Determining the fitting point from the run of the Brunt-Vaisala and Lamb frequencies
			 *****************************************************************************************/	
			for (k = 0; k<meshModel; k=k+2)
				if ((N2_model[k] > omega_r2_noDim) || (L2_model[k]*deg_l_dbl*(deg_l_dbl+1.0) < omega_r2_noDim))
				{
					csi_final = csiRel_model[k];		
					break;
				}					
			/*******************************************************************************
			 ***  Looping on the IMAGINARY Omega
			 *******************************************************************************/	
			for (loopOmega_i = 0; loopOmega_i < numStepsOmega_i; loopOmega_i++)
			{
				dbl_loopOmega_i = static_cast<double>(loopOmega_i);
				omega_i_noDim = (dbl_loopOmega_i*w_max_i_noDim + (numStepsOmega_i_dbl-dbl_loopOmega_i)*w_min_i_noDim)/numStepsOmega_i_dbl;
				
				if (detailedRiccati)
				{
					cout<< "#################################################################"<< endl;
					cout<< "########                                                 ########"<< endl;
					cout<< "######## calculating the evolution of the riccati matrix ########" << endl;
					cout<< "######## for the eigenfrequency:                         ########" << endl;
					if(tidesFlag)
					{
						cout<< "######## w tidal =  " << omega_r_noDim << "               ########" << endl;
						cout<< "######## w final =  " << csi_final << "               ########" << endl;
					}
					else
					{
						cout<< "######## w real  =  " << omega_r_noDim << "               ########" << endl;
						cout<< "######## w imag  =  " << omega_i_noDim << "               ########" << endl;
					}			
					cout<< "########                                                 ########"<< endl;
					cout<< "#################################################################"<< endl;;
				}
				
				/*Printing in a file the global properties if I am searching for an eigenfrequency*/
				if(detailedRiccati == 0)
					writeF_globalProperties(globalPropertiesSystem, M1_noDim*grams_units/Msun_toG, R1_noDim*cm_units/Rsun_toCm, L1_Lsun, Teff_K, 
											momInertia_noDim*grams_units*pow(cm_units,2.0)/(Msun_toG*pow(Rsun_toCm, 2.0)), M2_noDim*grams_units/Msun_toG, 
											Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, 
											a_noDim*cm_units/Rsun_toCm, ecc, omega_r_noDim, omega_i_noDim, PbreakUp_min, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
											epsilon_T_noDim, Clmk_noDim, Glmk_2_noDim, Glmk_3_noDim, klmk);
				/**************************************************************************
				 **************************************************************************
				 **************     Integrating the Riccati matrices in 2 steps:
				 **************     1) from surface to fitting point for n==0
				 **************     2) from center to fitting point for n==1
				 **************************************************************************
				 **************************************************************************/
				gsl_matrix_set_zero(matR_in_xFit);
				gsl_matrix_set_zero(matR_out_xFit);
				gsl_matrix_set_zero(matPermInfo_in);
				gsl_matrix_set_zero(matPermInfo_out);
				counterPermutIn = 0;
				counterPermutOut = 0;
				
				checkIn = 0; 
				checkOut = 0;
				for(n = 0; n < 2; n++)
				{
					/**************************************************************************************
					 ***  At the first time step I am integrating a permulation of R because of the BCs:
					 *** ==> create matrix T, and matTapplied to store all permutations that will be applied.
					 **************************************************************************************/
					pickInitialPermutationMatrix(matT, vecPermutInit_BC, rowsT, adiabatic, tidesFlag);
					gsl_matrix_memcpy(matTapplied, matT);
					/*************************************************************************************
					 ***   - Initial conditions for integrator (If adiabatic, Vi contain rows Riccati).
					 ***   - updating matPermInfo_#
					 ***   - indexing the permutations:	 the initial permutation applied has index 10000. 
					 *************************************************************************************/								
					if (n==0) //surface to fit.
					{
							InitialConditionsRij_surface(y_Riccati, rowsRic, adiabatic, deg_l_dbl, omega_r_noDim, omega_i_noDim, V_surf, 
														 delAD_surf, del_surf, epsilon_T_noDim, Clmk_noDim, tidesFlag, rho_surf_noDim, g_surf_noDim, 
														 surfaceRel);														
						
						csiRel = surfaceRel;
						h = -h_init;
						
						if(detailedRiccati)
						{
							for(k=0; k<sizeR; k++)
								gsl_vector_set(vec_y_Riccati_all,k, y_Riccati[k+position_Rprime_y_Riccati]);
							
							store_IC_ForIntegrationEigenfunctions(matPermInfo_in, csiRel, csi_final, 
																  vecPermutInit_BC, vec_y_Riccati_all, 
																  vecPermutInit_BC, counterPermutIn, sizeR);
							
							gsl_matrix_set(matPermInfo_in, 0, 2+rowsT+sizeR+rowsT, 10000);
						}//if(detailedRiccati)
					}//if (n==0)
					else //if (n==0)  //center to fit.
					{
						InitialConditionsRij_center(y_Riccati, rowsRic, adiabatic, deg_l_dbl, c1_center, omega_r_noDim, omega_i_noDim, tidesFlag);
						
						csiRel = centerRel;
						h = h_init;
						if(detailedRiccati)
						{
							
							for(k=0; k<sizeR; k++)
								gsl_vector_set(vec_y_Riccati_all,k, y_Riccati[k+position_Rprime_y_Riccati]);
							
							store_IC_ForIntegrationEigenfunctions(matPermInfo_out, csiRel, csi_final, 
																  vecPermutInit_BC, vec_y_Riccati_all, 
																  vecPermutInit_BC, counterPermutOut,sizeR);
							
							gsl_matrix_set(matPermInfo_out, 0, 2+rowsT+sizeR+rowsT, 10000);
							
						}//if(detailedRiccati)
					}//else								
					
					/* Calculating the norm of R */
					normR = 0.0;
					for (k=0; k<sizeR; k++) 
						normR = normR + pow(y_Riccati[k+position_Rprime_y_Riccati], 2.0);
					normR = sqrt(normR);					
					/*************************************************************************************
					 ***   Write on file the initial conditions, to make sure I have the starting point
					 ***   (in the while loop I don't print out the very first integration step)
					 *************************************************************************************/								
					if(nonAdiabatic)
					{
						if (n==0){OpenAndLabel_Rij_inOut("output/riccatiMatrices_in_", counterPermutIn, elmntsRin_f, rowsRic, tidesFlag);}
						else
							if(integratorFC){OpenAndLabel_Rij_inOut("output/riccatiMatrices_out_", counterPermutOut, elmntsRout_f, rowsRic, tidesFlag);}								
					}
					
					if(detailedRiccati)
					{					
						for(k=0; k<sizeR; k++)
							gsl_vector_set(vec_y_Riccati_all,k, y_Riccati[k+position_Rprime_y_Riccati]);
						
						if(n==0) //surface to fit
							writeF_Riccati_inOut(elmntsRin_f, csiRel, csiRel*cm_units/Rsun_toCm, vec_y_Riccati_all, 
												 sizeR, omega_r_noDim, omega_i_noDim, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
												 Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, tidesFlag);
						else  //center to fit
						{
							if(integratorFC)
								writeF_Riccati_inOut(elmntsRout_f, csiRel, csiRel*cm_units/Rsun_toCm, vec_y_Riccati_all, 
												 sizeR, omega_r_noDim, omega_i_noDim, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
												 Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, tidesFlag);
						}
						
					}//if(detailedRiccati)	
					/**************************************************************************
					 **************				Integrating
					 **************
					 ************** If eigenfunctions are calculated:
					 ************** - step1a) I store 3 consecutive interation info on Rij and the norm
					 ************** - step2a) I check wether the norm has a max. if it does I 
					 ************** set the riccati norm to the 1st elm of the 3 and re-start
					 ************** the integration from there, thus forcing a permutation when
					 ************** the norm has a peak.
					 ************** - step3a) After the new integration before permuting I check
					 ************** wether the new norm is bigger than the previous one. 
					 ************** If not, I re-integrate with a smaller timestep "h"
					 **************************************************************************/
					counterIntegrationSteps = 0;
					permutationApplied = 0;
					Nint_R = 0;
					Nint_R_tot = 0;
					Nint_R_atCheckNorm = 1;
					checkMaxNormR = 0;
					checkMaxNormR_OK = 1;
					
					if(detailedRiccati && nonAdiabatic && n==0){LimitToRiccati = 1.0e16;}
					
					if(n == 0)
						for(j = 0; j < 3; j++)
							for(i = 0; i < sizeR + 2 ; i++){gsl_matrix_set(matRelmnts_in,j, i, 0.0);}					
					
					if(adiabatic)
						sys_Riccati_R_toUse = sys_Riccati_R_adiabatic;
					else
						sys_Riccati_R_toUse = sys_Riccati_R_nonAdiabatic;										
					/*****************************************************************************************
					 ***** open file for riccati components and integrator state numbered based on the number of permutations
					 *****************************************************************************************/
					if(nonAdiabatic)
					{
						if(integratorRforV == 3)
						{
							if (n==0)							
								OpenAndLabel_rkf45_state("output/integratorState_in_", counterPermutIn, rkf45_state_in, rowsRic);
							else
								if(integratorFC)
									OpenAndLabel_rkf45_state("output/integratorState_out_", counterPermutOut, rkf45_state_out, rowsRic);							
						}						
					}

					/*****************************************************************************************
					 *****************************************************************************************/					
					while (fabs(csiRel - csi_final) > 0.0)
					{	  
						
						status = gsl_odeiv_evolve_apply  (e_Riccati, c_Riccati, s_Riccati, 
														  &sys_Riccati_R_toUse, &csiRel, csi_final, &h, y_Riccati);
						
						if (status != GSL_SUCCESS){break;}
												
						if(counterIntegrationSteps == 1)
						{			
							if(n == 0){surfaceFor_rkf45 = csiRel;}
							else{centerFor_rkf45 = csiRel;}
							if(integratorRforV == 3)
							{
								cout << "new star's boundaries = " << centerFor_rkf45 << "   " << surfaceFor_rkf45 << endl;
							}
						}//if(counterIntegrationSteps == 0)
						counterIntegrationSteps++;
						/************************************************************************
						 *** Extracting rkf45 state. This containts the coefficients of the fitting
						 *** functions describing the behavior of Rij between two consecutive integration
						 *** points,					 
						 *** 
						 *** - from center to fit I save everything
						 *** - from surface to fit, I save 3 consecutive elements following 
						 *** the way of saving and computing Rij described right above
						************************************************************************/
						if(integratorRforV == 3 && n == 1 && integratorFC)
							storeIntegratorState_out(vec_rkf45_state_out, s_Riccati, csiRel, sizeRij_calc);
						/************************************************************************/						
						Nint_R++;
						Nint_R_tot++;
						
						normR = 0.0;
						for (k = 0; k < sizeR; k++){normR = normR + pow(y_Riccati[k], 2.0);}
						normR = sqrt(normR);
						
						if(detailedRiccati && nonAdiabatic && n==0)
						{
							/*this is step (3a) described above*/		
							if(checkMaxNormR)
							{
								if(normR > normR_0){checkMaxNormR_OK = 1;}
								else
								{
									checkMaxNormR_OK = 0;
									Nint_R = Nint_R_atCheckNorm;
									
									/*re-setting initial conditions for integrator*/
									csiRel = gsl_matrix_get(matRelmnts_in, 0, 0);
									for (i = 0; i < sizeR; i++){y_Riccati[i] = gsl_matrix_get(matRelmnts_in, 0, 1 + i);}
									h = h/10.0;								
								}							
							}//if(checkMaxNormR)
							
							/*this is step 1a above*/
							if(checkMaxNormR_OK)
							{
								if(Nint_R <= 3)
								{
									/*waiting 3 steps to store at least the first 3 elements*/
									gsl_matrix_set(matRelmnts_in, Nint_R-1, 0, csiRel);					
									for(i = 1; i < sizeR + 1 ; i++){gsl_matrix_set(matRelmnts_in, Nint_R-1, i, y_Riccati[i-1]);}										
									gsl_matrix_set(matRelmnts_in, Nint_R-1, sizeR + 1, normR);
																		
									if(integratorRforV == 3)
										storeIntegratorState_in(mat_rkf45_state_in, s_Riccati, csiRel, sizeRij_calc, Nint_R-1);

								}//if(Nint_R <= 3)
								else
								{
									/*for each new element is move 1-->0 and 2-->1, and store 2.*/
									for(i = 0; i < sizeR + 2; i++)
									{
										gsl_matrix_set(matRelmnts_in, 0, i, gsl_matrix_get(matRelmnts_in, 1, i));
										gsl_matrix_set(matRelmnts_in, 1, i, gsl_matrix_get(matRelmnts_in, 2, i));
									}							
									gsl_matrix_set(matRelmnts_in, 2, 0, csiRel);					
									for(i = 1; i < sizeR + 1 ; i++){gsl_matrix_set(matRelmnts_in, 2, i, y_Riccati[i-1]);}										
									gsl_matrix_set(matRelmnts_in, 2, sizeR + 1, normR);
									
									if(integratorRforV == 3)
									{
										for(i = 0; i < size_rkf45_state; i++)
										{
											gsl_matrix_set(mat_rkf45_state_in, 0, i, gsl_matrix_get(mat_rkf45_state_in, 1, i));
											gsl_matrix_set(mat_rkf45_state_in, 1, i, gsl_matrix_get(mat_rkf45_state_in, 2, i));
										}//for(i = 0; i < size_rkf45_state; i++)
										storeIntegratorState_in(mat_rkf45_state_in, s_Riccati, csiRel, sizeRij_calc, 2);
									}//if(integratorRforV == 3)
									
								}//else
								
								if(Nint_R_atCheckNorm > 0 && (Nint_R == Nint_R_atCheckNorm + 1))
								{
									gsl_matrix_set(matRelmnts_in, 0, 0, csiRel);					
									for(i = 1; i < sizeR + 1 ; i++){gsl_matrix_set(matRelmnts_in, 0, i, y_Riccati[i-1]);}										
									gsl_matrix_set(matRelmnts_in, 0, sizeR + 1, normR);		
									
									
									if(integratorRforV == 3)
										storeIntegratorState_in(mat_rkf45_state_in, s_Riccati, csiRel, sizeRij_calc, 0);
								}
							}//if(checkMaxNormR_OK)	
							/*this is step 2a described above*/

							if(checkMaxNormR == 0 && Nint_R >= 3)
							{
								normR_0 = gsl_matrix_get(matRelmnts_in, 0, sizeR + 1);
								normR_1 = gsl_matrix_get(matRelmnts_in, 1, sizeR + 1);
								normR_2 = gsl_matrix_get(matRelmnts_in, 2, sizeR + 1);
								
								if((normR_1 > normR_0) && (normR_1 > normR_2))
								{
									LimitToRiccati = normR_0-0.1;
									checkMaxNormR = 1;
									Nint_R_atCheckNorm = Nint_R;
									
									/*re-setting initial conditions for integrator*/
									csiRel = gsl_matrix_get(matRelmnts_in, 0, 0);
									for (i = 0; i < sizeR; i++)
										y_Riccati[i] = gsl_matrix_get(matRelmnts_in, 0, 1 + i);	
									
									/*printing on file the initial integration step line*/
									for(i = 0; i < sizeR; i++)
										gsl_vector_set(vec_y_Riccati_all_test, i, gsl_matrix_get(matRelmnts_in, 0, i + 1));
									
									writeF_Riccati_inOut(elmntsRin_f, gsl_matrix_get(matRelmnts_in, 0, 0), 
														 gsl_matrix_get(matRelmnts_in, 0, 0)*cm_units/Rsun_toCm, vec_y_Riccati_all_test, 
														 sizeR, omega_r_noDim, omega_i_noDim, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
														 Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, tidesFlag);									
								
									if(integratorRforV == 3)
									{
										for (counterIntState = 0; counterIntState < size_rkf45_state; counterIntState++)
											gsl_vector_set(vec_rkf45_state_in, counterIntState, gsl_matrix_get(mat_rkf45_state_in, 0, counterIntState));
										
										writeF_integratorState_inOut(rkf45_state_in, vec_rkf45_state_in, size_rkf45_state);
									}
									
								
								}//if((normR_1 - normR_0 > 1e-12) && (normR_1 - normR_2 > 1.0e-12))								
							}//if(checkMaxNormR == 0 && Nint_R >= 3)
							
						}//if(detailedRiccati && nonAdiabatic && n==0)

						if(detailedRiccati && counterIntegrationSteps % 10000 == 0)
							cout << csiRel  << " " << normR << "  " << LimitToRiccati << ", #integrations = " << counterIntegrationSteps << "  h =  "<< h << "   " << Nint_R <<  "   " << Nint_R_atCheckNorm << "   " << checkMaxNormR_OK << endl; 
												
						/*To avoid numerical noise*/
						if(adiabatic)
							for(i = 0; i<sizeR+rowsRic; i++)
								if(fabs(y_Riccati[i]) > 0.0 && fabs(y_Riccati[i]) < 1.0e-16)
									y_Riccati[i] = 0.0;
						
						/* This is part of step (2a) described above. After the peak in the norm has been found and I restarted the integration
						 from one step before, I want to permute (Nint_R == Nint_R_atCheckNorm) is not satisfied*/
						if(detailedRiccati && nonAdiabatic && n == 0 && checkMaxNormR && Nint_R == Nint_R_atCheckNorm){;}
						else
						{
							if (permutationApplied)
							{
								if (n==0) {gsl_matrix_set(matPermInfo_in,counterPermutIn, 0, csiRel);}
								else {gsl_matrix_set(matPermInfo_out,counterPermutOut, 0, csiRel);}
								permutationApplied = 0;
							}											
							/*********************************
							 ***	Fill in R 
							 ********************************/											
							m = 0;
							for (k=0; k<rowsRic; k++)
								for (j=0; j<rowsRic; j++, m++)
									gsl_matrix_set(matR,k,j, y_Riccati[m + position_Rprime_y_Riccati]);		
							/***********************************
							 *** calculate ||R|| 
							 ***********************************/					
							normR = 0.0;
							for (k=0; k<rowsRic; k++)
								for (j=0; j<rowsRic; j++)
									normR = normR + pow(gsl_matrix_get(matR, k, j), 2.0);					
							normR = sqrt(normR);						
							/********************************************************
							 ***    If I'm calculating V
							 ***    - write on files the Rij (permuted and not)
							 ***    For the non-adiabatic case, since I am forcing a 
							 ***    permutation every time the norm has a maximum I have
							 ***    to modify the way I print out, to avoid having a non-monotonic r
							 ***    - update the matrix matPermInfo_#
							 *********************************************************/												
							if(detailedRiccati)
							{		
								/* store permuted elements R'# */													
								for(k=0; k<sizeR; k++)
									gsl_vector_set(vec_y_Riccati_all,k, y_Riccati[k+ position_Rprime_y_Riccati]);
								
								if(n==0) //surface to fit
								{									
									/**********************************************************************
									 ***    Write R' and update matPermInfo_#
									 **********************************************************************/												
									if(adiabatic)
										writeF_Riccati_inOut(elmntsRin_f, csiRel, csiRel*cm_units/Rsun_toCm, vec_y_Riccati_all, 
															 sizeR, omega_r_noDim, omega_i_noDim, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
															 Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, tidesFlag);
									else
									{
										if(Nint_R >=3)
										{
											for(i = 0; i < sizeR; i++)
												gsl_vector_set(vec_y_Riccati_all_test, i, gsl_matrix_get(matRelmnts_in, 0, i + 1));											
											
											writeF_Riccati_inOut(elmntsRin_f, gsl_matrix_get(matRelmnts_in, 0, 0), 
																 gsl_matrix_get(matRelmnts_in, 0, 0)*cm_units/Rsun_toCm, vec_y_Riccati_all_test, 
																 sizeR, omega_r_noDim, omega_i_noDim, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
																 Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, tidesFlag);
											
											if(integratorRforV == 3)
											{
												for (counterIntState = 0; counterIntState < size_rkf45_state; counterIntState++)
													gsl_vector_set(vec_rkf45_state_in, counterIntState, gsl_matrix_get(mat_rkf45_state_in, 0, counterIntState));
												
												writeF_integratorState_inOut(rkf45_state_in, vec_rkf45_state_in, size_rkf45_state);
											}//if(integratorRforV == 3)											
										}//if(Nint_R >=3)
									}//else
									
									savedRijIn++;									   
									
									update_csi_and_r_ij_forEndInterval(matPermInfo_in, csiRel, vec_y_Riccati_all, counterPermutIn, sizeR);
									
								}//if(n==0)
								else	//center to fit
								{
									if(integratorFC)
									{
										writeF_Riccati_inOut(elmntsRout_f, csiRel, csiRel*cm_units/Rsun_toCm, vec_y_Riccati_all, 
															 sizeR, omega_r_noDim, omega_i_noDim, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
															 Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, tidesFlag);
										
										if(integratorRforV == 3){writeF_integratorState_inOut(rkf45_state_out, vec_rkf45_state_out, size_rkf45_state);}
									}

									savedRijOut++;
									
									update_csi_and_r_ij_forEndInterval(matPermInfo_out, csiRel, vec_y_Riccati_all, counterPermutOut, sizeR);
									
								}//else
							}//if(detailedRiccati)	
							/***********************************************************************************************
							 *** 
							 ***  If ||R|| > user-defined limit:
							 ***  
							 ***  - Find the permutation with minimum ||R||, and store it in T.
							 ***  - Permute R, and re-set IC for integrator.
							 *** 
							 ***  If I'm calculating V:
							 *** 
							 ***  - Store in matPermInfo_# the info about T.
							 ***  - Multiply T with any other permutation applied before, and store it in matPermInfo_# 
							 ***  (this way we can always go back to the original R).
							 *** 
							 ***********************************************************************************************/						 
							if(normR > LimitToRiccati && (fabs(csiRel - csi_final) > 1.0e-16))
							{
								if(permutationNumber > MAX_NUMBER_PERMUTATIONS)
								{errorMessageAndExit("CAFein.c", "increase MAX_NUMBER_PERMUTATIONS");}
								
								if(detailedRiccati){cout << permutationNumber << " @ "  << csiRel << "-->cazzo permuting..." <<endl;}
								permutationApplied = 1;
								/*store in T the permutation with ||R_Real|| minimum. R_Real is unchanged from it.*/
								find_permutation_with_minNorm(matT, matT00, matT01, matT10, matT11, matR, matDummySizeR, matDummySizeR_2,
															  matDummySizeR_2, vecPermutIndices, vecDummySizeR2, 
															  vecDummySizeT, permDummySizeR, rowsRic, PermutationMinNormNumber, 
															  adiabatic, tidesFlag);
								
								/*recalculating R_real again, this time overwriting R*/
								permute_Riccati_matrix_R(matT, matR, matT00, matT01, matT10, matT11, 
														 matR, matDummySizeR, matDummySizeR_2, matDummySizeR_3, 
														 vecDummySizeR2, permDummySizeR, rowsRic);
																
								normRminimum = 0.0;
								for (k=0; k<rowsRic; k++) 
									for(j=0; j<rowsRic; j++)
										normRminimum = normRminimum + pow(gsl_matrix_get(matR,k,j), 2.0);
								normRminimum = sqrt(normRminimum);
								
								/* moving the Riccati limit if new permutation has bigger norm (not to keep permuting)*/
								if(normRminimum > LimitToRiccati){LimitToRiccati = normRminimum + normRminimum/10.0;}
								
								/*If i am computing the non-adiabatic eigenfunctions then I permute only when the norm reaches a maximum*/
								if(detailedRiccati && nonAdiabatic && n==0 && checkMaxNormR){LimitToRiccati = 1.0e16;}
								
								/*re-setting initial conditions for integrator*/
								m = 0;
								for (k=0; k<rowsRic; k++)
									for (j=0; j<rowsRic; j++, m++)
										y_Riccati[m+position_Rprime_y_Riccati] = gsl_matrix_get(matR, k, j);
								
								
								/*adding eigenfunctions part to the integrator if I am in the adiabatic case. 
								 At this stage I am integrating Riccati again (R00..R0N) 
								 (just not to have issues with the integrator if I set dV/dx = 0)*/
								if(adiabatic)
									for(k = 0; k<rowsRic; k++)
										y_Riccati[k+position_Vprime_y_Riccati] = y_Riccati[k+position_Rprime_y_Riccati];
								
								/*re-set integrator and  timestep h*/
								gsl_odeiv_evolve_reset(e_Riccati);
								gsl_odeiv_step_reset(s_Riccati);
								
								if(n==0){h = -h_init;}
								else {h = h_init;}		
																
								/*if I am calculating V*/
								if(detailedRiccati)
								{
									if(n==0){store_Tstep_ForIntegrationEigenfunctions(matPermInfo_in, vecPermutIndices, 
																					  sizeR, counterPermutIn, PermutationMinNormNumber);}																	
									else{store_Tstep_ForIntegrationEigenfunctions(matPermInfo_out, vecPermutIndices, 
																				  sizeR, counterPermutOut, PermutationMinNormNumber);}
								}//if(detailedRiccati)
								
								/*recalculating T considering all previous permutations (T = Tn Tn-1.......T1 T0 ) and updating 
								 the matrix storing all permutations applied */
								matrix_mul_rows_by_cols(matT, matTapplied, matT, rowsT, vecDummySizeT2);
								gsl_matrix_memcpy(matTapplied, matT);
								
								if(detailedRiccati)
								{	
									if(n==0){store_Ttot_ForIntegrationEigenfunctions(matT, rowsT, rowsT, matPermInfo_in, counterPermutIn);}
									else{store_Ttot_ForIntegrationEigenfunctions(matT, rowsT, rowsT, matPermInfo_out, counterPermutOut);}
									
								}//if(detailedRiccati)
								permutationNumber++;
								
								Nint_R = 0;
								Nint_R_atCheckNorm = 0;
								checkMaxNormR = 0;
								
								/***********************************************************************************************
								******* close the file containing the integrator state and open a new one
								 ***********************************************************************************************/								
								if(nonAdiabatic)
								{
									if (n==0)
									{
										elmntsRin_f.close();
										OpenAndLabel_Rij_inOut("output/riccatiMatrices_in_", counterPermutIn, elmntsRin_f, rowsRic, tidesFlag);
									}
									else
									{
										if(integratorFC)
										{
											elmntsRout_f.close();
											OpenAndLabel_Rij_inOut("output/riccatiMatrices_out_", counterPermutOut, elmntsRout_f, rowsRic, tidesFlag);
										}		
									}
									
								
									if(integratorRforV == 3)
									{
										if(n == 0)
										{
											rkf45_state_in.close();
											OpenAndLabel_rkf45_state("output/integratorState_in_", counterPermutIn, rkf45_state_in, rowsRic);									
										}
										else
											if(integratorFC)
											{
												rkf45_state_out.close();
												OpenAndLabel_rkf45_state("output/integratorState_out_", counterPermutOut, rkf45_state_out, rowsRic);
											}
									}
								}
								/***********************************************************************************************
								 ***********************************************************************************************/
							}//if (normR > LimitToRiccati)	
						}//else of if(nonAdiabatic && n == 0 && checkMaxNormR){;}
						/****************************************************************
						 ***  If I am at the fitting point:
						 ***  
						 ***  - Permute back to original R
						 ***  - store R in matR_in/out_xFit
						 *** 
						 *** matR_in and matR_out at xFit are needed to determine if the current
						 *** omega is an eigenfrequency, and to calculate the IC for V
						 ***  
						 ***  det(matR_out - matR_in)fit = 0 ==> omega is an eigenfrequency.
						 ****************************************************************/											
						if (fabs(csiRel - csi_final) <= 1.0e-16)
						{
							/*when almost at the surface just set it equal to the surface and end integration*/
							csiRel = csi_final;
							back_to_original_Riccati_R(matT, matT00, matT01, matT10, matT11, matR, matR, matDummySizeR, 
													   matDummySizeR_2, matDummySizeR_3, vecDummySizeR2, permDummySizeR, rowsRic);	
							
							if(n==0)//surface to fit
							{
								gsl_matrix_memcpy(matR_in_xFit, matR);
								checkIn++; 
							}
							else//center to fit
							{												
								gsl_matrix_memcpy(matR_out_xFit, matR);
								checkOut++; 
							}
						}//if (csi == csi_final)		
						/****************************************************************/
					}// while
					gsl_odeiv_evolve_reset(e_Riccati);
					gsl_odeiv_step_reset(s_Riccati);
										
					if(detailedRiccati && nonAdiabatic && n==0)
					{
						for(j = 1; j<=2; j++)
						{
							for(i = 0; i < sizeR; i++)
								gsl_vector_set(vec_y_Riccati_all_test,i, gsl_matrix_get(matRelmnts_in, j, i + 1));
							
							writeF_Riccati_inOut(elmntsRin_f, gsl_matrix_get(matRelmnts_in, j, 0), 
												 gsl_matrix_get(matRelmnts_in, j, 0)*cm_units/Rsun_toCm, vec_y_Riccati_all_test, 
												 sizeR, omega_r_noDim, omega_i_noDim, deg_l_dbl, deg_m_dbl, deg_k_dbl, 
												 Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, tidesFlag);
							
							if(integratorRforV == 3)
							{
								for (counterIntState = 0; counterIntState < size_rkf45_state; counterIntState++)
									gsl_vector_set(vec_rkf45_state_in, counterIntState, gsl_matrix_get(mat_rkf45_state_in, j, counterIntState));
								
								writeF_integratorState_inOut(rkf45_state_in, vec_rkf45_state_in, size_rkf45_state);
							}//if(integratorRforV == 3)
							
						}//for(j = 1; j<=2; j++)
					}//if(detailedRiccati && nonAdiabatic && n==0)
				
					if(nonAdiabatic)
					{
						if (n==0){elmntsRin_f.close();}
						else
							if(integratorFC){elmntsRout_f.close();}
						
						if(integratorRforV == 3)
						{
							if (n==0){rkf45_state_in.close();}							
							else
								if(integratorFC){rkf45_state_out.close();}
						}
					}
				}//for(n = nStart; n < nEnd; n++)
				
				if(detailedRiccati)
				{
					cout << "Intervals where permutations have been applied from surface to fit..." << endl;
					for(i=counterPermutIn ; i>0; i--)
						cout << gsl_matrix_get(matPermInfo_in,i,0) << "   " << gsl_matrix_get(matPermInfo_in,i,1) << endl;
				}
				
				
				if (checkIn != 1 || checkOut !=	1){errorMessageAndExit("CAFein.c", "csi_final not reached (checkIn/Out !=1)");}
				
				matrix_diff (matR_out_xFit, matR_in_xFit, matR_diff, rowsRic);
				
				det_RoMinusRin = determinantNbyN(matR_diff, rowsRic, permDummySizeR, matDummySizeR);
				
				/* vecDetDiff_r_i stores Re(determinant) and Im(determinant), which is only Re in the adiabatic case.*/
				gsl_vector_set(vecDetDiff_r_i,0,det_RoMinusRin);
				gsl_vector_set(vecDetDiff_r_i,1,0.0);				
				/************************************************************************************
				 ****** For R = 6x6 I have to extract the real and imaginary part of R to calculate 
				 ****** the condition on the determinant (recall the form of R from note #3 on top of
				 ****** this file)
				 ******
				 ****** Rreal = R00 = R11  = matR_diff_rr = matR_diff_ii
				 ****** Rimag = R10 = -R01 = matR_diff_ir = -matR_diff_ri
				 ************************************************************************************/
				if(nonAdiabatic)
				{
					extract_submatrices(matR_diff, matR_diff_rr, matR_diff_ri, matR_diff_ir, matR_diff_ii, static_cast<int>(rowsRic/2));					
					complexDeterminant3by3(matR_diff_rr, matR_diff_ir, vecDetDiff_r_i);								
				}
				
				writeF_Riccati_xFit(riccatiMcsi1, omega_r_noDim, omega_i_noDim, matR_in_xFit, matR_out_xFit, vecDetDiff_r_i, rowsRic, sizeDetVector);
				
				cout << "*** w_r = " << omega_r_noDim << " w_i = " << omega_i_noDim
				<< "  Re[det(Ro-Ri)] = " << gsl_vector_get(vecDetDiff_r_i,0) << " Im[det(Ro-Ri)] = " << gsl_vector_get(vecDetDiff_r_i,1) << endl;	
				
			}//for (loopOmega_i = 0; loopOmega_i < numStepsOmega_i; loopOmega_i++)
			
			/***************************************************************************
			 *** If secant method is applied, update the parameters for the calculation 
			 *** of the new omega_r_noDim and omega_i_noDim.
			 ***************************************************************************/			
			if(secant_omega_ri)
			{
				g_r_n_2  = g_r_n_1;
				g_i_n_2  = g_i_n_1;
				omegaR_n_2 = omegaR_n_1;
				omegaI_n_2 = omegaI_n_1;
				
				g_r_n_1  =  gsl_vector_get(vecDetDiff_r_i,0);
				g_i_n_1  =  gsl_vector_get(vecDetDiff_r_i,1);			
				omegaR_n_1 = omega_r_noDim;
				omegaI_n_1 = omega_i_noDim;
				
				cout << "delta(omegaR) =  " << fabs(omegaR_n_2  - omegaR_n_1) << endl;
				cout << "delta(omegaI) =  " << fabs(omegaI_n_2  - omegaI_n_1) << endl;
				
				/*************************************************************************************
				 *** keep applying the secant method till the difference between adjacent omega <1e-15
				 *** or both the real and imaginary determinant become close enough to zero
				 ************************************************************************************/							 
				if((fabs(omegaR_n_2 - omegaR_n_1) < 1.0e-15 && fabs(omegaI_n_2 - omegaI_n_1) < 1.0e-15) || 
				   ((fabs(gsl_vector_get(vecDetDiff_r_i,0)) > 0.0 && fabs(gsl_vector_get(vecDetDiff_r_i,0)) < 1.0e-15) &&				
					(fabs(gsl_vector_get(vecDetDiff_r_i,1)) > 0.0 && fabs(gsl_vector_get(vecDetDiff_r_i,1)) < 1.0e-15)))
				{
					cout << "....omega real..........omega imag..............det Real.........................det Imaginary " << endl;
					cout << omega_r_noDim <<	"   " << omega_i_noDim << "   " << gsl_vector_get(vecDetDiff_r_i,0) << "   " << gsl_vector_get(vecDetDiff_r_i,1) << endl;
				}
				else{loopOmega_r = -1;}
			}//if(secant_omega_ri)			
			counterIntStepsOmegaR_OmegaI0++;
			
			gsl_odeiv_evolve_reset(e_Riccati);
			gsl_odeiv_step_reset(s_Riccati);
		}//for (loopOmega_r=0; loopOmega_r<numStepsOmega_r; loopOmega_r++)		
		elmntsRin_f.close();
		if(integratorFC)
			elmntsRout_f.close();
		/**************************************************************************************************************************
		 **************************************************************************************************************************
		 **************************************************************************************************************************
		 **************************************************************************************************************************
		 **************
		 **************    Calculation of the eigenfunctions as described in section 2.3 of from Takata Loffler (2004)
		 **************    2004PASJ...56..645T
		 **************
		 **************    The configuration of U and V is as described in Note #1 above
		 **************
		 **************************************************************************************************************************
		 **************************************************************************************************************************
		 **************************************************************************************************************************
		 **************************************************************************************************************************/		
		if(detailedRiccati && (tidesFlag == 0))
		{
			cout<< "###################################################"<< endl;
			cout<< "########                                   ########"<< endl;
			cout<< "########         Rij calculated.           ########" << endl;
			cout<< "########     Finding the eigenfunction...  ########"<< endl;;
			cout<< "########                                   ########"<< endl;
			cout<< "###################################################"<< endl;;
		}
		if(tidesFlag)
		{
			cout<< "###################################################"<< endl;
			cout<< "########                                   ########"<< endl;
			cout<< "########         Rij calculated.           ########" << endl;
			cout<< "########       Finding the tidal           ########"<< endl;;
			cout<< "########          eigenfunction            ########"<< endl;;
			cout<< "########                                   ########"<< endl;
			cout<< "###################################################"<< endl;;
		}
		/******************************************************************************
		 *** Adiabatic: Reading Rij just calculated. These will be interpolated
		 *** with linear interpolation during the calculation of the eigenfunctions.
		 ******************************************************************************/
		if (adiabatic)
		{
			mesh_elements_Rin  = readInputFile(fileRiccatiInward.c_str(), elements_Rin);			
			if(integratorFC){mesh_elements_Rout = readInputFile(fileRiccatiOutward.c_str(), elements_Rout);}

			if (mesh_elements_Rin > mesh_elements_Rout){mesh_elements_Rmax = mesh_elements_Rin;}
			else{mesh_elements_Rmax = mesh_elements_Rout;}

			cout << "Mesh Rij Surface - fit ...= " << mesh_elements_Rin << endl;	
			cout << "Mesh Rij center - fit ....= " << mesh_elements_Rout << endl;	
			cout << "Max number of mesh points = " << mesh_elements_Rmax << endl;	
			
			if(mesh_elements_Rmax > MAX_NUMBER_INTEGRATIONSTEPS_AD)
				errorMessageAndExit("CAFein.c", "mesh_elements_Rmax > MAX_NUMBER_INTEGRATIONSTEPS_AD, increase MAX_NUMBER_INTEGRATIONSTEPS_AD!");		
		}//if (adiabatic)	

		counterPermut = 0;
		gsl_matrix_set_zero(matPermInfo);
		/******************************************************************************
		 *** allocating memory for the vector needed for the interpolation
		 *** of the Riccati elements and for steffen interpolant.
		 ******************************************************************************/		
		double **Rij_calc = NULL, *csi_calc = NULL, **fitRij_calcCoeffs = NULL, *Rij_calcFuncs = NULL,  
		**dR_calc_forSteffen = NULL, **h_forSteffen = NULL, **s_forSteffen = NULL, **rRel_and_rkf45_state = NULL, *Rij_from_rkf45_state = NULL; 

		if(detailedRiccati)
		{		
			/***********************************************
			 *re-initialize everything to zero just in case *
			 ***********************************************/
			for (i=0; i<4; i++)
			{
				rhoFuncs[i] = 0.0;
				gFuncs[i] = 0.0;
				VgFuncs[i] = 0.0;
				AstarFuncs[i] = 0.0;
				UFuncs[i] = 0.0;
				c1Funcs[i] = 0.0;
				VFuncs[i] = 0.0;
				delADFuncs[i] = 0.0;
				delFuncs[i] = 0.0;
				VtFuncs[i] = 0.0;
				ksFuncs[i] = 0.0;
				epsADFuncs[i] = 0.0;
				epsSFuncs[i] = 0.0;
				c2Funcs[i] = 0.0;
				c3Funcs[i] = 0.0;
				c4Funcs[i] = 0.0;
				dlnLR_dlnrFuncs[i] = 0.0;
				PFuncs[i] = 0.0;				
			}	
			/***********************************************************************************************
			 ***********************************************************************************************
			 **************       Integrating the Riccati matrices and eigenfunctions in 2 steps:
			 **************       1) from fitting point to surface for n==0
			 **************       2) from fitting point to centerRel for n==1
			 **************       The fitting point is where N^2 or L^2 == frequency considered and I 
			 **************       have no permutation to start with.
			 **************       
			 **************       The Rij just calculated can be either recalculated or read in input
			 **************       in the adiabatic case
			 ***********************************************************************************************
			 ***********************************************************************************************/			
			flagPrintOut = 0;
			string fileNameRiccatiPick, fileName_rkf45State_Pick;
			ifstream RiccatiElementsFile, rkf45_stateFile;

			for(n = nStart; n < nEnd; n++)
			{
				inwardOutward = n;
				/******************************************************************************
				 ******     - Picking the correct matrix matPermInfo_in (inward or outward).
				 ******     - Reading (R)_fit for BCs on eigefunctions.
				 ******     - Reading Rij calculated in case I want to interpolate
				 ******     - Reading the rkf45 state if rkf45 is used.
				 ******     
				 ******     In the second case:
				 ******		- if adiabatic: I am reading everything
				 ******		- if non adiabatic and no rkf45: I use simmetries in Note #3 on top of this file
				 ******		- if non adiabatic and rkf45: I read everything
				 ******     
				 ******     Note that:
				 ******     1) If rkf45 is used, the file containing k1 and k6 has ONE line
				 ******     less than the file containing Rij, because the star's boundaries are missing.
				 ******     				 
				 ******     2)When reading the Riccati elements and related parameters, I am reading
				 ******     them FROM the fitting point TO the star's boundaries.
				 ******     				 
				 ******     3)if rkf45 is used, here is how the numbering goes, given k computer below:
				 ******     (see also point (1) above)				 
				 ******     k ==  position of the fitting point in rkf45 state file
				 ******     k+1 ==  position of the fitting point in Rij file
				 ******     				 
				 ******     1 ==  position of the star's almost surface in rkf45 state file
				 ******     2 ==  position of the star's almost surface in Rij file
				 ******************************************************************************/				
				mesh_computedRij = (1-n)*mesh_elements_Rin + n*mesh_elements_Rout;
				
				counterPermut = counterPermutIn*(1-n) + n*counterPermutOut;
				
				if(adiabatic)
				{
					for (i=0; i<mesh_elements_Rmax; i++)
					{
						for (j = 0; j < sizeRij_calc; j++){Rij_calc_ad[i][j] = 0.0;}
						csi_calc_ad[i] = 0.0;						
					}

					/*reading riccati elements. The first two columns are the relative and absolute radius*/
					for (i = mesh_computedRij; i>=1; i--)
					{
						if(integratorFC)
						{
							csi_calc_ad[mesh_computedRij-i] = (1-n)*elements_Rin[(1-n)*i][0] + n*elements_Rout[n*i][0];							
							for (j = 0; j<sizeRij_calc; j++)
								Rij_calc_ad[mesh_computedRij-i][j] = (1-n)*elements_Rin[(1-n)*i][j + 2]+n*elements_Rout[n*i][j + 2];
						}
						else
						{
							csi_calc_ad[mesh_computedRij-i] = (1-n)*elements_Rin[(1-n)*i][0];							
							for (j = 0; j<sizeRij_calc; j++)
								Rij_calc_ad[mesh_computedRij-i][j] = (1-n)*elements_Rin[(1-n)*i][j + 2];							
						}
					}
				}//if(adiabatic)
				
				if (n==0) /* fit to surface */
				{
					gsl_matrix_memcpy(matPermInfo, matPermInfo_in);					
					gsl_matrix_memcpy(matRinit, matR_in_xFit);					
				}//if (n==0) 
				else /* fit to center */
				{
					gsl_matrix_memcpy(matPermInfo, matPermInfo_out);					
					gsl_matrix_memcpy(matRinit, matR_out_xFit);					
				}//else /* fit to center */
				/*************************************************************************************
				 ******     Looping on all the permutations applied during integration of R.
				 ******     For each permutation and integration interval I perform an integration
				 *************************************************************************************/
				for(z = counterPermut; z>=0; z--)
				{
					/*************************************************************************************
					 ******     non adiabatic case: opening files containing Riccati elements
					 ******     and integrator state.
					 *************************************************************************************/
					if(nonAdiabatic)
					{
						string permutationCntrLabel;
						ostringstream convert;								
						convert << z;
						permutationCntrLabel = convert.str();
						
						if (n==0)
							fileNameRiccatiPick = "output/riccatiMatrices_in_"+permutationCntrLabel+".dat";
						else
							fileNameRiccatiPick = "output/riccatiMatrices_out_"+permutationCntrLabel+".dat";

						mesh_elements_Rpick = countLinesInFile(fileNameRiccatiPick.c_str());															   
						cout << "Number of elements = " << mesh_elements_Rpick << endl;
						RiccatiElementsFile.open(fileNameRiccatiPick.c_str());
						
						if(integratorRforV == 3)
						{
							if (n==0)
								fileName_rkf45State_Pick = "output/integratorState_in_"+permutationCntrLabel+".dat";
							else
								fileName_rkf45State_Pick = "output/integratorState_out_"+permutationCntrLabel+".dat";

							rkf45_stateFile.open(fileName_rkf45State_Pick.c_str());
						}
					}//if(nonAdiabatic)
					
					/************************************************************************************
					 ***   Read the permutation Ttot and Tstep calculated before.
					 ***     
					 ***   (For the content of matPermInfo see explanation at the beginning of the code.
					 ************************************************************************************/
					for(i = 2; i<2+rowsT; i++)
						gsl_vector_set(vecPermutIndices,i-2, gsl_matrix_get(matPermInfo, z, i));
					
					fill_permutation_T(matTtot, vecPermutIndices, rowsT);	
					
					for(i = 2+rowsT+sizeR; i<(2+rowsT+sizeR+rowsT); i++)
						gsl_vector_set(vecPermutIndices,i-(2+rowsT+sizeR), gsl_matrix_get(matPermInfo, z, i));
					
					fill_permutation_T(matTstep, vecPermutIndices, rowsT);	
					/***************************************************************************
					 ***   The I.C on V at the fitting are given by Eq.(20) Takata Loffler 2004
					 ***   
					 ***   I then permute them based on Ttot computed during the calc. of Rij 
					 **********************************************************************/											
					if (z==counterPermut)
					{
						if((nonAdiabatic && tidesFlag) ||(adiabatic && tidesFlag == 0) )
						{
							/*finding the initial conditions on V via single value decomposition*/
							gsl_matrix_memcpy(matR_diff_SVD, matR_diff);
							dummy = gsl_linalg_SV_decomp(matR_diff_SVD, matV_SVD, vecS_SVD, work_SVD);
							
							/*The last column of matrix V(according to gsl manual, unhappy name...) contains the
							 initial conditions I want...*/
							for (i = 0; i<rowsRic; i++)
								cout << "SV" << i << " = " << gsl_vector_get(vecS_SVD, i) << endl;												
							
							for (i = 0; i<rowsRic; i++)
								gsl_vector_set(vecV, i, gsl_matrix_get(matV_SVD, i, rowsRic-1));												
						}
						else
							InitialConditionsUV_fittingPoint(matR_diff, vecV, adiabatic, tidesFlag);
						
						cout<< "###################################################"<< endl;
						cout<< "########                                   ########"<< endl;
						cout<< "######## Checking IC on V at fitting.....  ########" << endl;
						cout<< "########                                   ########"<< endl;
						cout<< "###################################################"<< endl;;
						
						for (int cazzo1=0; cazzo1< rowsRic; cazzo1++)
						{	
							checkIC_fittingPoint = 0.0;
							for (int cazzo2=0; cazzo2< rowsRic; cazzo2++)
							{
								checkIC_fittingPoint = checkIC_fittingPoint+
								gsl_matrix_get(matR_diff, cazzo1, cazzo2)*gsl_vector_get(vecV, cazzo2);
							}
							cout << "Rij*vi row" << cazzo1 << " = " << checkIC_fittingPoint << endl;							
						}
						
						/* calculate U */
						mul_matrix_by_vec(matRinit,vecV, vecU, vecU, rowsRic);
						
						for(i = 0; i < rowsRic; i++)
							cout << "v" << i << "_init = " <<  gsl_vector_get(vecV,i) << endl;						
						for(i = 0; i < rowsRic; i++)
							cout << "u" << i << "_init = " <<  gsl_vector_get(vecU,i) << endl;						

						/*Permute V with matTtot and over-write V, while U is still not permuted */
						permute_Riccati_V(vecU, vecV, matTtot, matT00, matT01, matT10, matT11, vecDummySizeR_1, vecDummySizeR_2, rowsRic);					
						
						/* Move V in matrix used for IC of integrator*/					
						gsl_vector_memcpy(vecVprime,vecV);
					}//if (z==counterPermut)	
					/**********************************************************************
					 ***    Read integration extrema.
					 ***   
					 ***    CsiIn_force and csiFin_force are passed to the integrator to force Rij 
					 ***    to be interpolated within the integration interval
					 **********************************************************************/					
					csiRel    = gsl_matrix_get(matPermInfo, z, 1);
					csi_final = gsl_matrix_get(matPermInfo, z, 0);	

					/*If rkf45 is used, the file containing rkf45 state has no boundaries, 
					 so I have to consider the last mesh point considered in rkf45 state*/
					if(nonAdiabatic)
					{
						if(integratorRforV == 3 && z == 0)
						{
							if(n == 0){csi_final = surfaceFor_rkf45;}
							else{csi_final = centerFor_rkf45;}
							cout << "new csi_final for rkf45 = " << csi_final << endl;
						}
					}
					double csi_surface = gsl_matrix_get(matPermInfo, 0, 0);

					csiIn_force  = csiRel; 
					csiFin_force = csi_final; 					
					/******************************************************************************
					 *** nonAdiabatic: Reading Rij just calculated. These will be interpolated
					 *** during the calculation of the eigenfunctions.
					 *** 
					 *** If Rij are computed using rkf45, then rkf45 status can be used
					 *** to determine the fitting function describing Rij across the whole integration
					 *** interval. 
					 *** In any other case, Rij SHOULD be interpolated using Steffen 
					 *** (apart for the adiabatic case where linear interpolation is good enough).		 
					 ******************************************************************************/					
					if(nonAdiabatic)
					{
						/*counting the number of computed Rij that are in the integration interval as given by the permutations applied*/
						mesh_elements_Rij = countElementsInInterval(RiccatiElementsFile,mesh_elements_Rpick, csiRel, csi_final, inwardOutward);
						cout << "permutation = " << z << ", Riccati elements in the radius interval = " << mesh_elements_Rij << endl;
						/*Allocate memory*/

						Rij_calc = dfunc(mesh_elements_Rij,sizeRij_calc);
						csi_calc = dfun(mesh_elements_Rij);	
						
						if(integratorRforV == 3)
						{
							rRel_and_rkf45_state = dfunc(mesh_elements_Rij, size_rkf45_state);
							Rij_from_rkf45_state = dfun(sizeRij_calc);
						}
						else
						{
							fitRij_calcCoeffs = dfunc(static_cast<int>(4*sizeRij_calc), mesh_elements_Rij);
							Rij_calcFuncs     = dfun(static_cast<int>(2*sizeRij_calc));
							dR_calc_forSteffen = dfunc(mesh_elements_Rij, sizeRij_calc); 
							h_forSteffen	   = dfunc(mesh_elements_Rij-1, sizeRij_calc); 
							s_forSteffen	   = dfunc(mesh_elements_Rij-1, sizeRij_calc); 			
						}	

						/*filling in the various vectors*/
						fillInRiccatiElements(RiccatiElementsFile, mesh_elements_Rpick, csiRel, csi_final, inwardOutward,
											  sizeRij_calc, csi_calc, Rij_calc, dummyTable);
						
						if(integratorRforV == 3)
							fillIn_rkf45_state(rkf45_stateFile, mesh_elements_Rpick, csiRel, csi_final, inwardOutward, 
											   size_rkf45_state, rRel_and_rkf45_state, dummyTable);
					
						/*now I can erase the file by opening it again and over-writing them with nothing*/
						RiccatiElementsFile.close();
						dummyFile.open(fileNameRiccatiPick.c_str());
						dummyFile.close();
						
						if(integratorRforV)
						{
							rkf45_stateFile.close();
							dummyFile.open(fileName_rkf45State_Pick.c_str());
							dummyFile.close();
						}
					}//if(nonAdiabatic)
					/************************************************************************************
					 ***  Structure containing the parameters for the integrator:
					 ***  OnwardOutward is any # >1 and is used only during the integration of the eigenfunctions.
					 **********************************************************************************/			
					params_integrator_V_nonAd_struct params_eigenfunc_nonAdiabatic = {csiIn_force, csiFin_force, rowsRic, meshModel, omega_r_noDim, 
						omega_i_noDim, deg_l_dbl, VgFuncs, AstarFuncs, UFuncs, c1Funcs, VFuncs, delADFuncs, delFuncs, VtFuncs, ksFuncs, 
						epsADFuncs, epsSFuncs, c2Funcs, c3Funcs, c4Funcs, dlnLR_dlnrFuncs, fitVgCoeffs, fitAstarCoeffs, fitUcoeffs, 
						fitC1coeffs, fitVCoeffs, fitDelADcoeffs, fitDelCoeffs, fitVtCoeffs, fitKsCoeffs, fitEpsADcoeffs, fitEpsScoeffs, 
						fitC2coeffs, fitC3coeffs, fitC4coeffs, fitdlnLR_dlnrCoeffs, csiRel_model, matA, matB, matC, matD, matCR, matR, matdR, 
						matTtot, matM, matDummySizeR, matDummySizeR_2, matDummySizeT, matDummySizeT_2, vecDummySizeR2, vecDummySizeT2, vecDummySizeR, 
						vecDummySizeR_3, permDummySizeT, inwardOutward, interpol_rij, sizeRij_calc, Rij_calc, csi_calc, fitRij_calcCoeffs, Rij_calcFuncs, 
						startInterpol_rij, endInterpol_rij, nonAdiabatic, tidesFlag};
					
					params_integrator_V_nonAd_rkf45_struct params_eigenfunc_nonAdiabatic_rkf45 = {csiIn_force, csiFin_force, rowsRic, meshModel, omega_r_noDim, 
						omega_i_noDim, deg_l_dbl, VgFuncs, AstarFuncs, UFuncs, c1Funcs, VFuncs, delADFuncs, delFuncs, VtFuncs, ksFuncs, 
						epsADFuncs, epsSFuncs, c2Funcs, c3Funcs, c4Funcs, dlnLR_dlnrFuncs, fitVgCoeffs, fitAstarCoeffs, fitUcoeffs, 
						fitC1coeffs, fitVCoeffs, fitDelADcoeffs, fitDelCoeffs, fitVtCoeffs, fitKsCoeffs, fitEpsADcoeffs, fitEpsScoeffs, 
						fitC2coeffs, fitC3coeffs, fitC4coeffs, fitdlnLR_dlnrCoeffs, csiRel_model, matA, matB, matC, matD, matCR, matR, matdR, 
						matTtot, matM, matDummySizeR, matDummySizeR_2, matDummySizeT, matDummySizeT_2, vecDummySizeR2, vecDummySizeT2, vecDummySizeR, 
						vecDummySizeR_3, permDummySizeT, inwardOutward, interpol_rij, sizeRij_calc, Rij_calc, csi_calc, 
						startInterpol_rij, endInterpol_rij, nonAdiabatic, tidesFlag, rRel_and_rkf45_state, Rij_from_rkf45_state, size_rkf45_state};
					
					gsl_odeiv_system sys_Riccati_R_eigenfunc_nonAdiabatic       = {func_Riccati_Vonly, jac_Riccati_Vonly, dimensionODE_V, &params_eigenfunc_nonAdiabatic};
					gsl_odeiv_system sys_Riccati_R_eigenfunc_nonAdiabatic_rkf45 = {func_Riccati_Vonly_rkf45, jac_Riccati_Vonly, dimensionODE_V, &params_eigenfunc_nonAdiabatic_rkf45};		
					/*********************************************************************
					 ***    Find the interval [csi_1, csi_2] where to interpolate R.
					 ***    I will pass this to the integrator to avoid looping across all
					 ***    the elements at each integration step.
					 ***    
					 ***    In the adiabatic case the integrator authomalically applies LI.					 
					 *********************************************************************/				
					endInterpol_rij = 0;
					startInterpol_rij = 0;				
					
					/*the case n==0 is for fit to surface*/
					/*the case n==1 is for fit to center*/
					if(adiabatic)
					{
						for(i=loopInterpol; i< (1-n)*mesh_elements_Rin + n*mesh_elements_Rout; i++)
						{
							if (((1-n)*csi_calc_ad[i] <= (1-n)*csiRel) && (n*csi_calc_ad[i] >= n*csiRel)){startInterpol_rij = i;}
							if (((1-n)*csi_calc_ad[i] <= (1-n)*csi_final) && (n*csi_calc_ad[i] >= n*csi_final)){endInterpol_rij = i;}
							else{break;}
						}
						loopInterpol = endInterpol_rij;
					}
					else
					{
						startInterpol_rij = 0;
						endInterpol_rij = mesh_elements_Rij-1;
						
						if (integratorRforV == 3 && n == 0 && z == 0)
						{
							for (i = 0; i<mesh_elements_Rij; i++)
							{
								if (csi_calc[i] > csi_final)
								{
									endInterpol_rij = i-1;
									cout << "I just found = " << endInterpol_rij << endl;
									break;
								}
							}
						}
						cout << "where it ends now = " << endInterpol_rij << endl;

					}

					if (tidesFlag && z == 0 && fabs(csiRel-csi_final) < 1.0e-5)
					{
						cout << "\n******************************************" << endl;
						cout << "**********                          ******" << endl;
						cout << "**********  Skipping permutation:   ******" << endl;
						cout << "**********           "<<counterPermut<<"      ******"<< endl;
						cout << "**********     r_f - r_i < 1e-5     ******" << endl;
						cout << "********** permutation at surface   ******" << endl;
						cout << "******************************************\n" << endl;										
						break;
					}							
					
					/*calculating Steffen interpolants. This will be used in the non-adiabatic case only.
					 In the adiabatic case the integrator authomatically applies linear interpolation*/
					if(nonAdiabatic && integratorRforV !=3)
					{
						for (i = 0; i < static_cast<int>(4*sizeRij_calc); i++)
							for (j = 0; j < mesh_elements_Rij; j++)
								fitRij_calcCoeffs[i][j] = 0.0;

						calcSteffenInterp_RiccatiElmnts(sizeRij_calc, startInterpol_rij, endInterpol_rij, dR_calc_forSteffen, 
														h_forSteffen, s_forSteffen, csi_calc, Rij_calc, fitRij_calcCoeffs);
					}
					/*************************************************
					 ***    Initializing integrator variables 
					 ***    - adiabatic case:
					 ***    If R is interpolated, then the structure is
					 ***    
					 ***    y[0] = v0, y[1] = v1, 
					 ***    y[2] = v0, y[3] = v1, 
					 ***    y[4] = v0, y[5] = v1, 
					 ***    
					 ***    If R is integrated, then the structure is
					 ***    
					 ***    y[0] = r00, y[1] = r01, 
					 ***    y[2] = r10, y[3] = r11, 
					 ***    y[4] = v0,  y[5] = v1, 
					 ***    
					 ***    - non adiabatic case:
					 ***    
					 ***    y[0] = v0, y[1] = v1, y[2] = v2
					 ***    y[3] = v3, y[4] = v4, y[5] = v5
					 ***    
					 *************************************************/
					/*Rij as read from the stored permutation info */
					if(adiabatic)
						for (i = 0; i < sizeR; i++)
							y_Riccati_V[i+position_Rprime_y_Riccati] = gsl_matrix_get(matPermInfo,z, i+position_Rprime_matPermutations) ;					
					
					/*Vi*/
					for (i = 0; i < rowsRic; i++)
						y_Riccati_V[i + position_Vprime_y_Riccati] = gsl_vector_get(vecVprime, i);										
					
					/*if Rij are interpolated, I fill in the remaining variables with Vi */
					if(interpol_rij && adiabatic)
						for (i = 0; i <rowsRic; i++)
							for (j = 0; j <rowsRic; j++)
								y_Riccati_V[rowsRic*i+j+position_Rprime_y_Riccati] = y_Riccati_V[j+position_Vprime_y_Riccati];
					
					gsl_odeiv_evolve_reset(e_eigenfunc);
					gsl_odeiv_step_reset(s_eigenfunc);					
					/***********************************
					 ***      Print on screen stuff
					 ***********************************/				
					if(tidesFlag == 0)
					{
						if(n==0 && flagPrintOut==0)
						{
							cout << "\n******************************************" << endl;
							cout << "**********" << endl;
							cout << "**********     xFit to surface:" << endl;
							cout << "**********   " << counterPermut <<  " permutations" << endl;
							cout << "******************************************\n" << endl;					
							flagPrintOut++;
						}
						if(n==1 && flagPrintOut==1)
						{
							cout << "\n******************************************" << endl;
							cout << "**********" << endl;
							cout << "**********     xFit to center:" << endl;
							cout << "**********   " << counterPermut <<  " permutations" << endl;
							cout << "******************************************\n" << endl;					
							flagPrintOut--;
						}
						cout << "\nPermutation # " << counterPermut - z << "\n"<< endl;	
						cout << "**************************************" << endl;
						cout << "***** Integration extremes: **********" << endl;
						cout << "**************************************" << endl;
						cout << "*** x_init = " << csiRel << endl;
						cout << "*** x_fin .= " << csi_final << endl;
						cout << "**************************************" << endl;
					}//if(tidesFlag == 0)
					/**************************************************************************
					 **************************************************************************
					 **************				Integrating
					 **************		
					 **************	 I am calculating V' (and R' if it is not interpolated). 
					 **************	 From V', R' I calculate U' = R'V'		 
					 **************		
					 **************************************************************************
					 **************************************************************************/									
					IntegrationStepsEigenfunc = 0;
					
					/*if rkf45 is used, Steffen is authomatically not used*/
						
					if (n==0){h = h_init_eigenfunc;}
					else{h = -h_init_eigenfunc;}					
					
					/* consider two different integrator for adiabatic or non adiabatic equations. In the adiabatic case I have the option
					 of integrating R and V simultaneously*/
					if(adiabatic){sys_Riccati_R_toUse = sys_Riccati_R_eigenfunc_adiabatic;} 
					else
					{
						if (integratorRforV == 3)
							sys_Riccati_R_toUse = sys_Riccati_R_eigenfunc_nonAdiabatic_rkf45;
						else
							sys_Riccati_R_toUse = sys_Riccati_R_eigenfunc_nonAdiabatic;
					}
					
					//FIXME to erase if running with WD:
					if (n==0 && nonAdiabatic){csi_final = csi_calc[endInterpol_rij-1];}
                    
                    
					while (fabs(csiRel - csi_final)>0)
					{	 		
						if (inwardOutward < 0 || inwardOutward > 1)
							errorMessageAndExit("CAFein.c", "inwardOutward MUST be 0 or 1 when integrating eigenfunctions");												
						
						status = gsl_odeiv_evolve_apply  (e_eigenfunc, c_eigenfunc, s_eigenfunc, &sys_Riccati_R_toUse, 
														  &csiRel, csi_final, &h, y_Riccati_V);
						
						if (status != GSL_SUCCESS){break;}
						
						IntegrationStepsEigenfunc++;
						
						
						/*To avoid numerical noise*/
						if(adiabatic)
						{
							for(i = 0; i<rowsRic; i++)
								if(fabs(y_Riccati_V[i]) > 0.0 && fabs(y_Riccati_V[i]) < 1.0e-16){y_Riccati_V[i] = 0.0;}
							
							/* If R' is calculated I read the corresponding y_Riccati_V and store them in the matrix matDummySizeR */
							/* R' is calculated only in the adiabatic case.*/
							m = 0;
							for (i = 0; i < rowsRic; i++)
								for (j = 0; j < rowsRic; j++, m++)
									gsl_matrix_set(matRprime, i, j, y_Riccati_V[m+position_Rprime_y_Riccati]);
						}//if(adiabatic)
                        
						/* if R' is interpolated I read matrix R which was filled in the integrator
						 R' can be interpolated in the adiabatic case and for sure in the non adiabatic case*/
						if (interpol_rij){gsl_matrix_memcpy(matRprime,matR);}
												
						/* read V'*/
						for (i = 0; i < rowsRic; i++)
							gsl_vector_set(vecVprime, i, y_Riccati_V[i + position_Vprime_y_Riccati]);
						
						/* calculating U' = R'V'*/
						mul_matrix_by_vec(matRprime, vecVprime, vecUprime, vecUprime, rowsRic);
						
						/* copy U' and V' in U and V, for permutation purposes..*/
						gsl_vector_memcpy(vecU, vecUprime);
						gsl_vector_memcpy(vecV, vecVprime);
						
						/*for output purposes, go back to original U and V overwriting vecUoriginal, and vecVoriginal*/
						back_to_original_Riccati_UV(matTtot, matT00, matT01, matT10, matT11, vecU, vecV, 
													vecUoriginal, vecVoriginal, matDummySizeT, matDummySizeT_2, 
													vecDummySizeR_1, vecDummySizeR_2, vecDummySizeT, permDummySizeT, 
													rowsRic);
						
                        
						/* For stability, the equations to integrate have been modified to integrate y_i/(x^(l-2)).
						 For output purposes I now go back to y_i*/
						for (i = 0; i < rowsRic; i++)
						{
							gsl_vector_set(vecVoriginal, i, gsl_vector_get(vecVoriginal,i) * pow(csiRel, deg_l_dbl-2.0));
							gsl_vector_set(vecUoriginal, i, gsl_vector_get(vecUoriginal,i) * pow(csiRel, deg_l_dbl-2.0));
						}
						
						evalSteffenInterp(meshModel,csiRel_model,fitRhoCoeffs,csiRel, rhoFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fit_g_Coeffs,csiRel, gFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fitPcoeffs,csiRel, PFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fitGamma1coeffs,csiRel, Gamma1Funcs);	
						evalSteffenInterp(meshModel,csiRel_model,fitN2coeffs,csiRel, N2Funcs);	
						evalSteffenInterp(meshModel,csiRel_model,fitTcoeffs,csiRel, TFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fit_Cp_Coeffs,csiRel, CpFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fitScoeffs,csiRel, SFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fit_Lr_Coeffs,csiRel, LrFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fitDelADcoeffs,csiRel, delADFuncs);	
						evalSteffenInterp(meshModel,csiRel_model,fitVCoeffs,csiRel, VFuncs);	
						
						if(numStepsPspin == 1 && ((IntegrationStepsEigenfunc < 100) || (IntegrationStepsEigenfunc % writeEigenfunctionsSkip == 0)))
						{
							
							/*Uncomment this one if you want all the Riccati coefficients together with the eigenfunctions*/
//							columnsEigenfuncFile = writeF_eigenfunctions(eigenfunctions, csiRel, csiRel*R1_Rsun, omega_r_noDim, omega_i_noDim, deg_l_dbl, 
//																		 rhoFuncs[0], gFuncs[0], vecUoriginal, vecVoriginal, matRprime, 
//																		 rowsRic, adiabatic, tidesFlag, Pspin_minutes, Porb_minutes, deg_m_dbl, deg_k_dbl,
//																		 tidesTimescale_a, GRtimescale_a, tidesTimescale_Omega, arg_Flmk*180.0/PI, Mod_Flmk);
							
							columnsEigenfuncFile = writeF_eigenfunctions_and_unperturbed_model(eigenfunctions, csiRel, csiRel*cm_units/Rsun_toCm, 						  																			
																							   vecUoriginal, vecVoriginal, rhoFuncs[0], gFuncs[0], 
																							   PFuncs[0], Gamma1Funcs[0], N2Funcs[0], delADFuncs[0], 
																							   TFuncs[0], CpFuncs[0], SFuncs[0], LrFuncs[0], 
																							   tidesTimescale_a, tidesTimescale_e, tidesTimescale_Omega,	   
																							   GRtimescale_a, arg_Flmk*180.0/PI, Mod_Flmk, 
																							   Pspin_noDim*sec_units/secPerMin, rowsRic, adiabatic, tidesFlag, linearityViol, 
																							   Jorb_noDim, Jspin_noDim, JdotOrb_noDim, JdotSpin_noDim, VFuncs[0]);
						}
						if ((IntegrationStepsEigenfunc % 1000 == 0) && numStepsPspin == 1 && nonAdiabatic)
							cout << "star's radius = " << csiRel << "....csi_final = " << csi_final  << "   " <<csi_calc[endInterpol_rij] << "   " <<gsl_vector_get(vecVoriginal, 0) << endl ;
						
						/*because the tidal eigenfunctions might have a point where they drop to zero close to the surface for numerical
						 reasons, I stop to 1e-5 from the surface.*/
						if(n == 0 && fabs(csiRel-csi_surface)<1.0e-5){csiRel = csi_final;}							
						
					}// while (end integration)

					/* permute U and V for calculating the initial conditions for the next integration step using matTstep
					 this time I overwrite U and V, although it doesn't matter at this point */
					back_to_original_Riccati_UV(matTstep, matT00, matT01, matT10, matT11, vecU, vecV, vecU, vecV, matDummySizeT, matDummySizeT_2, 
												vecDummySizeR_1, vecDummySizeR_2, vecDummySizeT, permDummySizeT, rowsRic);
					
					/* Re-set the initial conditions for V' */
					gsl_vector_memcpy(vecVprime, vecV);
					gsl_odeiv_evolve_reset(e_eigenfunc);
					gsl_odeiv_step_reset(s_eigenfunc);
					if (n==0){h = h_init_eigenfunc;}
					else{h = -h_init_eigenfunc;}	
					
					if(nonAdiabatic)
					{
						free(csi_calc);	
						for (i = 0; i < mesh_elements_Rij; i++){free(Rij_calc[i]);}
						free(Rij_calc);	
						
						if(integratorRforV == 3)
						{
							for (i = 0; i < mesh_elements_Rij; i++){free(rRel_and_rkf45_state[i]);}
							free(rRel_and_rkf45_state);	
							free(Rij_from_rkf45_state);
						}
						else
						{
							for(i = 0; i < 4*sizeRij_calc; i++){free(fitRij_calcCoeffs[i]);}
							free(fitRij_calcCoeffs);
							free(Rij_calcFuncs);		
							
							for (i = 0; i < mesh_elements_Rij; i++){free(dR_calc_forSteffen[i]);}
							free(dR_calc_forSteffen);
							
							for(i = 0; i < 	mesh_elements_Rij-1; i++)
							{
								free(h_forSteffen[i]); 
								free(s_forSteffen[i]); 	
							}
							free(h_forSteffen); 
							free(s_forSteffen); 	
						}	
						
					}//if(nonAdiabatic)										
				}//for(z = counterPermutIn; z>=0; z--)

				loopInterpol = 0;
				
				/*************************************************************************************
				 *************************************************************************************
				 ********		
				 ********	Normalize the eigenfunctions at the surface if tides are ON. 
				 ********	Now I have a file containing my tidal eigenfunctions not normalized.
				 ********	I want to normalize them so that y8 = 1 at the star's surface.
				 ********	Specifically: y8_real = 1 and y8_imag = 0.
				 ********	(see routine normalizeEigenfunctions_y8to1atR for an explanation)
				 ********		
				 *************************************************************************************
				 *************************************************************************************/									
				
				if(tidesFlag && n == 0)
				{					
					/*Extracting y8 at the surface*/
					
					if(adiabatic)
					{
						y8r_atR = gsl_vector_get(vecVoriginal, rowsRic-1);
						y8i_atR = 0.0;
					}
					else
					{
						y8r_atR = gsl_vector_get(vecVoriginal, (rowsRic/2)-1);
						y8i_atR = gsl_vector_get(vecVoriginal, rowsRic-1);
					}					
					GSL_SET_COMPLEX(&y_complex, y8r_atR, y8i_atR);
					
					arg_y8_atR     = gsl_complex_arg(y_complex);
					modulus_y8_atR = gsl_complex_abs(y_complex);
					
					/*************************************************************************************
					 *************************************************************************************
					 ********		
					 ********	Looping on the eigenfunctions horizontally to normalize them 
					 ********	AT THE SURFACE. I will write them in a separate file.
					 ********		
					 *************************************************************************************
					 *************************************************************************************/									
					for(i = 0; i< adiabatic*rowsRic + nonAdiabatic*rowsRic/2; i++)
					{
						/*extract the real and imaginary part of the eigenfunctions for U and V*/
						yReU = gsl_vector_get(vecUoriginal, i);
						yReV = gsl_vector_get(vecVoriginal, i);						
						yImU = 0.0;
						yImV = 0.0;
						
						if(nonAdiabatic)
						{
							yImU = gsl_vector_get(vecUoriginal, i + (rowsRic/2));
							yImV = gsl_vector_get(vecVoriginal, i + (rowsRic/2));						
						}
						
						/*extract modulus and argument of each eigenfunction for U*/
						GSL_SET_COMPLEX(&y_complex, yReU, yImU);
						arg_y_atR     = gsl_complex_arg(y_complex);
						modulus_y_atR = gsl_complex_abs(y_complex);
						
						cout << "arg_y_atR y8 = " << arg_y8_atR << ", arg_y_atR = " <<  arg_y_atR << "   ";
						/*normalizing U*/
						yReNormU = (modulus_y_atR/modulus_y8_atR)*(cos(arg_y_atR)*cos(arg_y8_atR)+sin(arg_y_atR)*sin(arg_y8_atR));
						yImNormU = (modulus_y_atR/modulus_y8_atR)*(sin(arg_y_atR)*cos(arg_y8_atR)-sin(arg_y8_atR)*cos(arg_y_atR));
						
						/*extract modulus and argument of each eigenfunction for V*/
						GSL_SET_COMPLEX(&y_complex, yReV, yImV);
						arg_y_atR     = gsl_complex_arg(y_complex);
						modulus_y_atR = gsl_complex_abs(y_complex);
						
						cout << "arg_y_atR = " <<  arg_y_atR << endl;
						
						/*normalizing V*/
						yReNormV = (modulus_y_atR/modulus_y8_atR)*(cos(arg_y_atR)*cos(arg_y8_atR)+sin(arg_y_atR)*sin(arg_y8_atR));
						yImNormV = (modulus_y_atR/modulus_y8_atR)*(sin(arg_y_atR)*cos(arg_y8_atR)-sin(arg_y8_atR)*cos(arg_y_atR));
						
						/*sobstitute the normalized eigenfunctions in the vector vecU/Voriginal_normalized*/
						gsl_vector_set(vecUoriginal_normalized, i, yReNormU);
						gsl_vector_set(vecVoriginal_normalized, i, yReNormV);
						
						if(nonAdiabatic)
						{
							gsl_vector_set(vecUoriginal_normalized, i + (rowsRic/2), yImNormU);
							gsl_vector_set(vecVoriginal_normalized, i + (rowsRic/2), yImNormV);
						}
					}//for(i = 0; i<(rowsRic/2) + (adiabatic*rowsRic/2); i++)
										
					/*************************************************************************************
					 *************************************************************************************
					 ********		
					 ********	To compute the timescales I need the total perturbation of the potential
					 ********		
					 *************************************************************************************
					 *************************************************************************************/									
					
					/* extracting y3, which contains the total perturbation of the potential y3 = Psi_T/(gr)*/
					gsl_vector_set(Psi_T_realImag, 0, gsl_vector_get(vecVoriginal_normalized, 0));
					gsl_vector_set(Psi_T_realImag, 1, 0);
					if(nonAdiabatic){gsl_vector_set(Psi_T_realImag, 1, gsl_vector_get(vecVoriginal_normalized, rowsRic/2));}
					
					evalSteffenInterp(meshModel,csiRel_model,fit_g_Coeffs,csiRel, gFuncs);	
					calculate_Flmk(deg_l_int, deg_m_int, deg_k_int, ecc, R1_noDim, a_noDim, epsilon_T_noDim, gFuncs[0], csiRel, Flmk_realImag, Psi_T_realImag);
					
					GSL_SET_COMPLEX(&Flmk,gsl_vector_get(Flmk_realImag, 0), gsl_vector_get(Flmk_realImag, 1));
					
					arg_Flmk = gsl_complex_arg(Flmk);
					Mod_Flmk = gsl_complex_abs(Flmk);
					
					/*if the argument is negative, move it between 0 and 2PI by adding 2PI
					 whenever it is negative*/
					if(arg_Flmk < 0.0){arg_Flmk = arg_Flmk + 2.0*PI;}
					/*in the adiabatic case since everything is real assume sinGamma in eq. (54) of Bart's paper to be 1*/
					if(adiabatic)
					{
						arg_Flmk = PI/2.0;
						cout << "adiabatic sin (gamma) = " << sin(arg_Flmk) << endl;
					}
					
					/*************************************************************************************
					 ******** Computing the dimensionless derivatives
					 *************************************************************************************/														
					dadt_tides_noDim = calculate_dAdT_tides(deg_l_int, deg_m_int, deg_k_int, ecc, R1_noDim, a_noDim, epsilon_T_noDim, Porb_noDim, 
															M1_noDim, M2_noDim, adiabatic, arg_Flmk, Mod_Flmk);					
					
					c_light_noDim = c_cgs/(cm_units/sec_units);
					dadt_GR_noDim = calculate_dAdT_GR(a_noDim, M1_noDim, M2_noDim, c_light_noDim);					
					
					dedt_tides_noDim = calculate_dedT_tides(deg_l_int, deg_m_int, deg_k_int, ecc, R1_noDim, a_noDim, epsilon_T_noDim, Porb_noDim, 
															M1_noDim, M2_noDim, adiabatic, arg_Flmk, Mod_Flmk);
					
					dOmegadt_tides_noDim = calculate_dOmegadT_tides(deg_l_int, deg_m_int, deg_k_int, ecc, R1_noDim, a_noDim, epsilon_T_noDim, Porb_noDim, 
																	M1_noDim, M2_noDim, adiabatic, momInertia_noDim, arg_Flmk, Mod_Flmk);
					/*************************************************************************************
					 ******** Computing the timescales and putting the dimensions back
					 *************************************************************************************/														
					/* Orbital separation */
					tidesTimescale_a = a_noDim/dadt_tides_noDim;
					tidesTimescale_a = tidesTimescale_a*sec_units/secPerYears;
					
					/* eccentricity */
					if(fabs(ecc - 0.0) <= 1.0e-10){tidesTimescale_e = 0.0;}
					else 
					{
						tidesTimescale_e = ecc/dedt_tides_noDim;
						tidesTimescale_e = tidesTimescale_e*sec_units/secPerYears;
					}
					
					/* star's spin */
					tidesTimescale_Omega = omegaSpin_noDim/dOmegadt_tides_noDim;
					tidesTimescale_Omega = tidesTimescale_Omega*sec_units/secPerYears;
					
					/* GR */
					GRtimescale_a  = a_noDim/dadt_GR_noDim;
					GRtimescale_a = GRtimescale_a*sec_units/secPerYears;
					
					/*************************************************************************************
					 ******** Computing the dimensionless angular momenta
					 *************************************************************************************/														
					Jorb_noDim = sqrt((M1_noDim+M2_noDim))*(M1_noDim*M2_noDim/(M1_noDim+M2_noDim))*sqrt(a_noDim*(1.0-ecc*ecc));
					
					Jspin_noDim = momInertia_noDim * omegaSpin_noDim;
					
					JdotOrb_noDim  = Jorb_noDim * 0.5 * ((dadt_tides_noDim/a_noDim) - (2.0*ecc/(1.0 - ecc*ecc))*dedt_tides_noDim);
					JdotSpin_noDim = Jspin_noDim * dOmegadt_tides_noDim/omegaSpin_noDim;
					
					cout << "JtotDot/Jtot = " << (JdotOrb_noDim+JdotSpin_noDim)/(Jorb_noDim+Jspin_noDim) << "   " << JdotOrb_noDim/Jorb_noDim << "   " << JdotSpin_noDim/Jspin_noDim << endl;
					
					/*
					 radial component of the tidal displacement with the inverse of the radial wavenumber*/
					/*******************************************************************************************************
					 *** check a possible violation of linearity by comparing the amplitude of the oscillatory part of the 
					 *******************************************************************************************************/																			
					kh2 = deg_l_dbl*(deg_l_dbl+1.0)/(csiRel*csiRel);
					Cs2 = Gamma1Funcs[0]*PFuncs[0]/rhoFuncs[0];
					kr2 = (N2Funcs[0] - omega_r_noDim*omega_r_noDim)*(Cs2*kh2 - omega_r_noDim*omega_r_noDim)/(omega_r_noDim*omega_r_noDim*Cs2);
					
					if(tidesFlag)
					{
						/* xi_r_dyn = y1*csiRel , xi_r_st = -y3*csiRel*/
						xi_r_dyn_R = gsl_vector_get(vecUoriginal_normalized, 0)*csiRel;
						xi_r_st_R = -gsl_vector_get(vecVoriginal_normalized, 0)*csiRel;
						
						if(nonAdiabatic)
						{
							xi_r_dyn_I = gsl_vector_get(vecUoriginal_normalized, rowsRic/2)*csiRel;
							xi_r_st_I = -gsl_vector_get(vecVoriginal_normalized, rowsRic/2)*csiRel;
						}
						modXi_r_osc = sqrt((xi_r_dyn_R - xi_r_st_R)*(xi_r_dyn_R - xi_r_st_R)+(xi_r_dyn_I - xi_r_st_I)*(xi_r_dyn_I - xi_r_st_I));
						if(modXi_r_osc > (1.0/sqrt(kr2))){linearityViol = 1;}
					}
					
					/*Uncomment this one if you want all the Riccati coefficients together with the eigenfunctions*/					
//					columnsEigenfuncFile = writeF_eigenfunctions(eigenfunctionsSurface, csiRel, csiRel*R1_Rsun, omega_r_noDim, omega_i_noDim, deg_l_dbl, 
//																 rhoFuncs[0], gFuncs[0], vecUoriginal_normalized, vecVoriginal_normalized, matRprime, 
//																 rowsRic, adiabatic, tidesFlag, Pspin_minutes, Porb_minutes, deg_m_dbl, deg_k_dbl, 
//																 tidesTimescale_a, GRtimescale_a, tidesTimescale_Omega, arg_Flmk*180.0/PI, Mod_Flmk);
					
					if (adiabatic && tidesFlag){Mod_Flmk = gsl_vector_get(Flmk_realImag, 0);}
					columnsEigenfuncFile = writeF_eigenfunctions_and_unperturbed_model(eigenfunctionsSurface, csiRel, csiRel*cm_units/Rsun_toCm, 						  																			
																					   vecUoriginal_normalized, vecVoriginal_normalized, rhoFuncs[0], gFuncs[0], 
																					   PFuncs[0], Gamma1Funcs[0], N2Funcs[0], delADFuncs[0], 
																					   TFuncs[0], CpFuncs[0], SFuncs[0], LrFuncs[0], 
																					   tidesTimescale_a, tidesTimescale_e, tidesTimescale_Omega,	   
																					   GRtimescale_a, arg_Flmk*180.0/PI, Mod_Flmk, 
																					   Pspin_noDim*sec_units/secPerMin, rowsRic, adiabatic, tidesFlag, linearityViol, 
																					   Jorb_noDim, Jspin_noDim, JdotOrb_noDim, JdotSpin_noDim, VFuncs[0]);
				}//if(tidesFlag && n == 0)

			}//for(n = nStart; n < nEnd; n++)
			
			/*Write the global properties of the system in a file*/
			writeF_globalProperties(globalPropertiesSystem, M1_noDim*grams_units/Msun_toG, R1_noDim*cm_units/Rsun_toCm, L1_Lsun, Teff_K, 
									momInertia_noDim*grams_units*pow(cm_units, 2.0)/(Msun_toG*pow(Rsun_toCm, 2.0)), M2_noDim*grams_units/Msun_toG, 
									Pspin_noDim*sec_units/secPerMin, Porb_noDim*sec_units/secPerMin, a_noDim*cm_units/Rsun_toCm, ecc, omega_r_noDim, 
									omega_i_noDim, PbreakUp_min, deg_l_dbl, deg_m_dbl, deg_k_dbl, epsilon_T_noDim, Clmk_noDim, Glmk_2_noDim, Glmk_3_noDim, klmk);
			
		}//if(detailedRiccati)

		/*******************************************
		 ***		cleaning up the memory
		 *******************************************/
		elements_Rin.clear();
		elements_Rout.clear();
		rRel_and_rkf45_stateIn.clear();
		rRel_and_rkf45_stateOut.clear();
		
		/*******************************************/
		/*******************************************/
	}//for(loopPspin = 0; loopPspin < numStepsPspin; numStepsPspin++)
	eigenfunctions.close();
	eigenfunctionsSurface.close();
	riccatiMcsi1.close();
	
	/**************************************************************************
	 **************************************************************************
	 **** 
	 **** If tides are used and only one loop on Pspin is applied, then
	 **** normalize the eigenfunctions so that y8 = 1 at the surface 
	 **** 
	 **************************************************************************
	 **************************************************************************/				
	if(tidesFlag && numStepsPspin == 1)
	{
		cout << "columnsEigenfuncFile = " << columnsEigenfuncFile << endl;
		
		/* Opening the new file that will contain the normalized eigenfunctions */
		const char *fileNameEigenfunctions_y81atSurf;
		string fileEigenfunctions_y81atSurf;
		
		
		fileEigenfunctions_y81atSurf = "output/eigenfunctions_y81atSurf.dat";
		
		fileNameEigenfunctions_y81atSurf = fileEigenfunctions_y81atSurf.c_str();
		
		normalizeEigenfunctions_y8to1atR(fileEigenfunctions.c_str(), fileNameEigenfunctions_y81atSurf, adiabatic,
										 rowsRic, tidesFlag, arg_y8_atR, modulus_y8_atR, columnsEigenfuncFile);		
		
	}//if(tidesFlag && numStepsPspin == 1)
	
	
	globalPropertiesSystem.close();
	
	free(csi_calc_ad);	
	
	for (i = 0; i < MAX_NUMBER_INTEGRATIONSTEPS_AD; i++)
		free(Rij_calc_ad[i]);
	free(Rij_calc_ad);	
	
	free(rhoFuncs);
	free(gFuncs);
	free(VgFuncs);
	free(AstarFuncs);
	free(UFuncs);
	free(c1Funcs);
	free(VtFuncs);
	free(VFuncs);
	free(delFuncs);
	free(delADFuncs);
	free(ksFuncs);
	free(c2Funcs);
	free(c3Funcs);
	free(c4Funcs);
	free(epsADFuncs);
	free(dlnLR_dlnrFuncs);
	free(epsSFuncs);
	free(PFuncs);
	free(Gamma1Funcs);
	free(N2Funcs);
	free(TFuncs);
	free(CpFuncs);
	free(SFuncs);
	free(LrFuncs);
	
	free(rho_model);
	free(g_model);
	free(csiRel_model);
	free(mass_model);
	free(Vg_model);
	free(Astar_model);
	free(U_model);
	free(c1_model);
	free(N2_model);
	free(L2_model);
	free(Vt_model);
	free(V_model);
	free(del_model);
	free(delad_model);
	free(ks_model);
	free(c2_model);
	free(c3_model);
	free(c4_model);
	free(epsAd_model); 
	free(dlnLR_dlnr_model);
	free(epsS_model);
	free(deltaRadRegion_model);
	free(P_model);
	free(Gamma1_model);
	free(T_model);
	free(Cp_model);
	free(S_model);
	free(Lr_model);
	
	
	for (i=0; i<4; i++)
	{
		free(fitRhoCoeffs[i]);
		free(fit_g_Coeffs[i]);
		free(fitVgCoeffs[i]);
		free(fitAstarCoeffs[i]);
		free(fitUcoeffs[i]);
		free(fitC1coeffs[i]);
		free(fitVtCoeffs[i]);
		free(fitVCoeffs[i]);
		free(fitDelCoeffs[i]);
		free(fitDelADcoeffs[i]);
		free(fitKsCoeffs[i]);
		free(fitC2coeffs[i]);
		free(fitC3coeffs[i]);
		free(fitC4coeffs[i]);
		free(fitEpsADcoeffs[i]);
		free(fitdlnLR_dlnrCoeffs[i]);
		free(fitEpsScoeffs[i]);
		free(fitPcoeffs[i]);
		free(fitGamma1coeffs[i]);
		free(fitTcoeffs[i]);
		free(fit_Cp_Coeffs[i]);
		free(fitScoeffs[i]);
		free(fit_Lr_Coeffs[i]);
		free(fitN2coeffs[i]);
	}
	
	free(fitRhoCoeffs);
	free(fit_g_Coeffs);
	free(fitVgCoeffs);
	free(fitAstarCoeffs);
	free(fitUcoeffs);
	free(fitC1coeffs);
	free(fitVtCoeffs);
	free(fitVCoeffs);
	free(fitDelCoeffs);
	free(fitDelADcoeffs);
	free(fitKsCoeffs);
	free(fitC2coeffs);
	free(fitC3coeffs);
	free(fitC4coeffs);
	free(fitEpsADcoeffs);
	free(fitdlnLR_dlnrCoeffs);
	free(fitEpsScoeffs);
	free(fitPcoeffs);
	free(fitGamma1coeffs);
	free(fitTcoeffs);
	free(fit_Cp_Coeffs);
	free(fitScoeffs);
	free(fit_Lr_Coeffs);
	free(fitN2coeffs);
	free(y_Riccati_all);
	
	
	gsl_odeiv_evolve_free (e_Riccati);
	gsl_odeiv_control_free (c_Riccati);
	gsl_odeiv_step_free (s_Riccati);
	gsl_odeiv_evolve_free (e_eigenfunc);
	gsl_odeiv_control_free (c_eigenfunc);
	gsl_odeiv_step_free (s_eigenfunc);
	
	gsl_matrix_free(matDummySizeT);
	gsl_matrix_free(matDummySizeT_2);
	gsl_matrix_free(matDummySizeR);
	gsl_matrix_free(matDummySizeR_2);
	gsl_matrix_free(matDummySizeR_3);
	gsl_matrix_free(matRoriginal);
	gsl_matrix_free(matT);
	gsl_matrix_free(matTapplied);
	gsl_matrix_free(matM);
	gsl_matrix_free(matT00);
	gsl_matrix_free(matT01);
	gsl_matrix_free(matT10);
	gsl_matrix_free(matT11);
	gsl_matrix_free(matR);
	gsl_matrix_free(matdR);
	gsl_matrix_free(matA);
	gsl_matrix_free(matB);
	gsl_matrix_free(matC);
	gsl_matrix_free(matD);
	gsl_matrix_free(matCR);
	gsl_matrix_free(matR_in_xFit);
	gsl_matrix_free(matR_out_xFit);
	gsl_matrix_free(matPermInfo_in);
	gsl_matrix_free(matPermInfo_out);
	gsl_matrix_free(matRelmnts_in);
	gsl_vector_free(vecPermutIndices);
	gsl_vector_free(vec_y_Riccati_all);
	gsl_vector_free(vec_y_Riccati_all_test);
	gsl_vector_free(vecPermutInit_BC);
	gsl_vector_free(vecDummySizeT);
	gsl_vector_free(vecDummySizeR2);
	gsl_vector_free(vecDummySizeR);
	gsl_vector_free(vecDummySizeR_3);
	gsl_vector_free(vecDummySizeT2);
	gsl_vector_free(vecDetDiff_r_i);
	gsl_permutation_free(permDummySizeR);
	gsl_permutation_free(permDummySizeT);
	gsl_matrix_free(matR_diff);
	gsl_matrix_free(matR_diff_SVD);
	gsl_matrix_free(matV_SVD);
	gsl_matrix_free(matR_diff_rr);
	gsl_matrix_free(matR_diff_ii);
	gsl_matrix_free(matR_diff_ri);
	gsl_matrix_free(matR_diff_ir);
	gsl_vector_free(vecV);
	gsl_vector_free(vecU);
	gsl_vector_free(vecVoriginal);
	gsl_vector_free(vecUoriginal);
	gsl_matrix_free(matRinit);
	gsl_matrix_free(matRprime);
	gsl_vector_free(vecVprime);
	gsl_vector_free(vecUprime);
	gsl_matrix_free(matTstep);
	gsl_matrix_free(matTtot);
	gsl_vector_free(vecDummySizeR_1);
	gsl_vector_free(vecDummySizeR_2);
	gsl_vector_free(vecVoriginal_normalized);
	gsl_vector_free(vecUoriginal_normalized);
	gsl_vector_free(Flmk_realImag);
	gsl_vector_free(Psi_T_realImag);
	gsl_vector_free(vecS_SVD);
	gsl_vector_free(work_SVD);

	if(integratorRforV ==3)
	{
		gsl_vector_free(vec_rkf45_state_in);
		gsl_vector_free(vec_rkf45_state_out);
		gsl_matrix_free(mat_rkf45_state_in);
	}//if(integratorRforV ==3)
	

	cout << "cazzo DONE" << endl;
	
	return 0;
}

