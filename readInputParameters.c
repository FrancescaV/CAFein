#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ostream>
#include <iomanip>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "readInputParameters.h"
#include "IOfiles.h"

using namespace std;

/*************************************************************************************
 ************ 
 ************	The routine described below read the input parameters from acc_and_dim.dat
 ************ 
 **************************************************************************************/

void readInputParams(const char *fileName, 
					 int &adiabatic_func, int &detailedRiccati_func,
					 int &secant_omega_ri_func, int &numStepsOmega_r_func, int &numStepsOmega_i_func, 
					 int &interpol_rij_func, int &reset_h_eachNsteps_func,
					 int &keep_hConst_forNsteps_func, int &writeEigenfunctionsSkip_func, int &integratorRforV_func, int &integratorV_func, 
					 int &integratorFC_func, int &integratorFS_func,
					 int &nStart_func, int &nEnd_func, double &l_func, double &h_init_func, double &eps_abs_func, 
					 double &eps_rel_func, double &LimitToRiccati_func, double &w_min_r_func, double &w_max_r_func, 
					 double &displaceEfreq_func, double &w_min_i_func, double &w_max_i_func, double &eigenv_r_func, 
					 double &eigenv_i_func, double &h_init_eigenfunc_func, double &eps_abs_eigenfunc_func, 
					 double &eps_rel_eigenfunc_func,
					 double &omegaR_n_2_func, double &omegaR_n_1_func, double &omegaI_n_2_func, double &omegaI_n_1_func, 
					 double &g_r_n_2_func, double &g_r_n_1_func, double &g_i_n_2_func, double &g_i_n_1_func, double &g_r_func, 
					 double &g_i_func, int &tidesFlag_func, int &numStepsPspin_func, double &deg_m_func, double &deg_k_func, 
					 double &M2_Msun_func, double &Porb_func, double &Pspin_min_func, double &Pspin_max_func, 
					 int &deg_l_int_func, int &deg_m_int_func, int &deg_k_int_func, double &ecc_func)
{
	typedef vector<double> Row;
	vector<Row> table;

	double w_min_r_tempo = 0.0, w_max_r_tempo = 0.0, 
	w_min_i_tempo = 0.0, w_max_i_tempo = 0.0, 
	Pspin_max_tempo = 0.0, Pspin_min_tempo = 0.0;
	
	int numStepsOmega_r_tempo = 0, numStepsOmega_i_tempo = 0, numStepsPspin_tempo = 0;
	
	
	int dummy = 0, addLine = 0;
	dummy = readInputFile(fileName, table);

	adiabatic_func = static_cast<int>(table[0 + addLine][1]);
	if(adiabatic_func > 1){errorMessageAndExit("readInputParameters.c", "adiabatic or not?!!");}

	//reading element 2, or #3 in acc_and_dim.dat
	
	addLine = addLine+3;	
	l_func		  		= static_cast<double>(table[1 + addLine][1]);
	deg_l_int_func      = static_cast<int>(table[1 + addLine][1]);
	h_init_func			= table[2 + addLine][1];
	eps_abs_func		= table[3 + addLine][1];
	eps_rel_func		= table[4 + addLine][1];
	LimitToRiccati_func	= table[5 + addLine][1];
	
	addLine = addLine+5;	
	detailedRiccati_func = static_cast<int>(table[6 + addLine][1]),
	secant_omega_ri_func = static_cast<int>(table[7 + addLine][1]);
	
	if(detailedRiccati_func> 1){errorMessageAndExit("readInputParameters.c", "unkown detailedRiccati!");}	
	if(secant_omega_ri_func > 1){errorMessageAndExit("readInputParameters.c", "unkown secant_omega_ri!");}	
	if(secant_omega_ri_func){detailedRiccati_func = 0;}
	
	if(adiabatic_func && secant_omega_ri_func)
		errorMessageAndExit ("readInputParameters.c", "secant_omega_ri ==1,adiabatic ==1!Secant method not applicable in adiabatic case.");
	
	/* Read the extremes for looping on omega_r and omega_i*/
	addLine = addLine+3;	
	w_min_r_tempo		  = table[8 + addLine][1];
	w_max_r_tempo		  = table[9 + addLine][1];
	numStepsOmega_r_tempo = static_cast<int>(table[10 + addLine][1]);
	displaceEfreq_func	  = table[11 + addLine][1];

	addLine = addLine+3;	
	w_min_i_tempo		  = table[12 + addLine][1];
	w_max_i_tempo		  = table[13 + addLine][1];
	numStepsOmega_i_tempo = static_cast<int>(table[14 + addLine][1]);
		
	addLine = addLine+5;	
	eigenv_r_func			= table[15 + addLine][1];
	eigenv_i_func           = table[16 + addLine][1];
	h_init_eigenfunc_func	= table[17 + addLine][1];
	eps_abs_eigenfunc_func	= table[18 + addLine][1];
	eps_rel_eigenfunc_func	= table[19 + addLine][1];
	interpol_rij_func		= static_cast<int>(table[20 + addLine][1]);

	addLine = addLine+5;
	if(secant_omega_ri_func)
	{		
		omegaR_n_2_func = table[21 + addLine][1];
		omegaR_n_1_func	= table[22 + addLine][1];
		omegaI_n_2_func	= table[23 + addLine][1];
		omegaI_n_1_func	= table[24 + addLine][1];
		
		g_r_n_2_func = table[25 + addLine][1];
		g_r_n_1_func = table[26 + addLine][1];
		g_i_n_2_func = table[27 + addLine][1];
		g_i_n_1_func = table[28 + addLine][1];
		
		g_r_func = table[29 + addLine][1];
		g_i_func = table[30 + addLine][1];
		
	}	
	
	addLine = addLine+5;
	//reading element 39, or #40 in acc_and_dim.dat
	writeEigenfunctionsSkip_func = static_cast<int>(table[31 + addLine][1]);
 	
	addLine = addLine+8;	
	integratorRforV_func		= static_cast<int>(table[32 + addLine][1]);
	integratorV_func			= static_cast<int>(table[33 + addLine][1]);
	
	/*Forcing bsimp for non-adiabatic eigenfrequencies*/
	if(adiabatic_func == 0 && detailedRiccati_func == 0)
	{
		integratorRforV_func = 2;	
	}
	if(adiabatic_func)
	{
		integratorRforV_func = 1;
		integratorV_func = 1;
	}
	integratorFC_func			= static_cast<int>(table[34 + addLine][1]);
	integratorFS_func			= static_cast<int>(table[35 + addLine][1]);
	
	if(integratorRforV_func > 3|| integratorRforV_func < 0) 
		errorMessageAndExit("readInputParameters.c", "unkown integratorRforV!");
	
	if(integratorV_func > 3 || integratorV_func < 0)
		errorMessageAndExit("readInputParameters.c", "unkown integratorV!");
	
	if(integratorFC_func > 1)
		errorMessageAndExit("readInputParameters.c", "unkown integratorFC!");
	
	if(integratorFS_func > 1)
		errorMessageAndExit("readInputParameters.c", "unkown integratorFS!");
	
	if(integratorFC_func == 0  && integratorFS_func == 1)
	{
		nStart_func = 0;
		nEnd_func = 1;
	}
	else if(integratorFC_func == 1  && integratorFS_func == 0)
	{
		nStart_func = 1;
		nEnd_func = 2;
	}
	else if(integratorFC_func == 1  && integratorFS_func == 1)
	{
		nStart_func = 0;
		nEnd_func = 2;
	}
	else {errorMessageAndExit("readInputParameters.c", "unkown integation interval!");}
	
	addLine = addLine+5;	
	tidesFlag_func			= static_cast<int>(table[37 + addLine][1]);

	Pspin_min_tempo			= table[40 + addLine][1];
	Pspin_max_tempo			= table[41 + addLine][1];
	numStepsPspin_tempo	    = static_cast<int>(table[42 + addLine][1]);

	if(tidesFlag_func)
	{
		M2_Msun_func   = table[38 + addLine][1];
		Porb_func	   = table[39 + addLine][1];
		deg_m_func	   = static_cast<double>(table[43 + addLine][1]);
		deg_m_int_func = static_cast<int>(table[43 + addLine][1]);

		deg_k_func		= static_cast<double>(table[44 + addLine][1]);
		deg_k_int_func  = static_cast<int>(table[44 + addLine][1]);
		ecc_func        = table[45 + addLine][1];
	}//if(tidesFlag_func)
	if(tidesFlag_func > 1)
		errorMessageAndExit("readInputParameters.c", "unkown tidesFlag!");	
	
	/***********************************************************
	 ***** 
	 ***** Handling the extremes of the loops on the various
	 ***** frequencies.
	 ***** 
	 ***** THe options are:
	 ***** 
	 ***** adiabatic...........Tides.......detailed riccati
	 *****		1				0			0........only calculating the natural eigenfrequencies
	 *****		1				0			1........only calculating the natural eigenfunctions
	 *****		1				1			1........only calculating the forced eigenfunctions
	 *****		0				0			0........only calculating the natural eigenfrequencies
	 *****		0				0			1........only calculating the natural eigenfunctions
	 *****		0				1			1........only calculating the forced eigenfunctions
	 ***********************************************************/

	if(adiabatic_func)
	{
		/*in the adiabatic case omega is only real*/
		w_min_i_func = 0.0;
		w_max_i_func = 0.0;
		numStepsOmega_i_func = 1;
		
		if(tidesFlag_func == 0)
		{
			/* without tides skip the loop on the spin */
			Pspin_max_func = 0.0;
			Pspin_min_func = 0.0;
			numStepsPspin_func = 1;
						
			if(detailedRiccati_func == 0)
			{
				/* loop on omega_r*/
				w_min_r_func = w_min_r_tempo;
				w_max_r_func = w_max_r_tempo;
				numStepsOmega_r_func = numStepsOmega_r_tempo;
			}//if(detailedRiccati_func == 0)
			else
			{
				/*Read eigenfrequencies in input and kill the loop on omega_r*/
				w_min_r_func = eigenv_r_func;
				w_max_r_func = eigenv_r_func;
				numStepsOmega_r_func = 1;				
			}//else detailed Riccati			
		}//if(tidesFlag_func == 0)
		else
		{
			/* with tides set up the loop on the spin */
			detailedRiccati_func = 1;
			Pspin_max_func = Pspin_max_tempo;
			Pspin_min_func = Pspin_min_tempo;
			numStepsPspin_func = numStepsPspin_tempo;
			numStepsOmega_r_func = 1;				
			/*and the loop on omega_r is calculated in the code*/
		}//else adiabatic with tides
	}//if(adiabatic_func)
	else
	{
		if(tidesFlag_func == 0)
		{
			/* without tides skip the loop on the spin */
			Pspin_max_func = 0.0;
			Pspin_min_func = 0.0;
			numStepsPspin_func = 1;			

			if(detailedRiccati_func == 0)
			{
				/* loop on omega_r and omega_i*/
				w_min_r_func = w_min_r_tempo;
				w_max_r_func = w_max_r_tempo;
				numStepsOmega_r_func = numStepsOmega_r_tempo;
				w_min_i_func = w_min_i_tempo;
				w_max_i_func = w_max_i_tempo;
				numStepsOmega_i_func = numStepsOmega_i_tempo;				
			}//if(detailedRiccati_func == 0)
			else
			{
				/*Read eigenfrequencies in input and kill the loop on omega_r and omega_i*/
				w_min_r_func = eigenv_r_func;
				w_max_r_func = eigenv_r_func;
				numStepsOmega_r_func = 1;
				w_min_i_func = eigenv_i_func;
				w_max_i_func = eigenv_i_func;
				numStepsOmega_i_func = 1;				
			}//else detailed Riccati			
		}//if(tidesFlag_func == 0)
		else
		{
			/* with tides set up the loop on the spin and kill the loop on omega_i*/
			detailedRiccati_func = 1;
			Pspin_max_func = Pspin_max_tempo;
			Pspin_min_func = Pspin_min_tempo;
			numStepsPspin_func = numStepsPspin_tempo;
			
			w_min_i_func = 0.0;
			w_max_i_func = 0.0;
			numStepsOmega_i_func = 1;
			numStepsOmega_r_func = 1;				
			/*and the loop on omega_r is calculated in the code*/
			
		}//else non adiabatic with tides
	}//else non adiabatic
	
	
	/* During the integration of V, Rij are interpolated by default in the non adiabatic case*/
	if(adiabatic_func == 0){interpol_rij_func = 1;}
	
	table.clear();
	
}
