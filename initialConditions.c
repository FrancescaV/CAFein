#include <iostream>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include "initialConditions.h"
#include "RiccatiMatrices_operations.h"

using namespace std;

const double PI = 2.0*acos(0.0);
/*************************************************************************************
 ************ 
 ************	The routines described below handle the initial conditions .
 ************   for eigenfrequencies (riccati components) and eigenfunctions
 **************************************************************************************/


void InitialConditionsRij_surface(double *y_Riccati_func, int rowsR_func, int adiabatic_func, double l_func, 
								  double omega_r_func, double omega_i_func, 
								  double V_func, double delAD_func, double del_func, double epsilon_t_func, 
								  double Clmk_func, int tidesFlag_func, double rho_surf_func, double g_surf_func, 
								  double surfaceRel_func)
{	
	int k = 0;
	double
	wR2 = omega_r_func*omega_r_func, 
	wI2 = omega_i_func*omega_i_func,
	wIwR = omega_i_func*omega_r_func,
	lPlus1 = 1.0 + l_func,
	coeff2 = (l_func*lPlus1)/((wI2 + wR2)*(wI2 + wR2)), 	
	coeff3 = -wI2 + wR2, 	
	coeff4 = (-4.0 + (-1.0 + coeff2)*coeff3 + V_func)*(-4.0 + (-1.0 + coeff2)*coeff3 + V_func) + 
	4.0*(1.0 + coeff2)*(1.0 + coeff2)*wI2*wR2, 
	delAD_V = delAD_func*V_func, 
	twoMinus4DelAdV = 2.0 - 4.0*delAD_V,
	z1 = (-4.0 + (l_func*lPlus1/wR2) - wR2)/V_func,
	z2 = (-1.0 - l_func + (l_func*lPlus1/wR2))/V_func, 
	onePlusZ1 = 1.0 + z1, 
	onePlusZ2 = 1.0 + z2;

	
	for(k = 0; k < rowsR_func*rowsR_func; k++){y_Riccati_func[k] = 0;}
	
	if(adiabatic_func)
	{
		if(tidesFlag_func)
		{
			y_Riccati_func[0] = 1.0;													//r00	
			y_Riccati_func[1] = -1.0;													//r01
			y_Riccati_func[2] = 0.0;													//r02
			
			y_Riccati_func[3] = -4.0*PI*rho_surf_func/g_surf_func;						//r10	
			y_Riccati_func[4] = -l_func -1.0 + (4.0*PI*rho_surf_func/g_surf_func);		//r11	
			y_Riccati_func[5] = -epsilon_t_func*(2.0*l_func + 1.0)*Clmk_func/g_surf_func;//r12
			
			y_Riccati_func[6] = 0.0;													//r20
			y_Riccati_func[7] = 0.0;													//r21
			y_Riccati_func[8] = 1.0;													//r22

			y_Riccati_func[9]  = y_Riccati_func[0];										//v0
			y_Riccati_func[10] = y_Riccati_func[1];										//v1
			y_Riccati_func[11] = y_Riccati_func[2];										//v2
		}//if(tidesFlag_func)
		else
		{
			/* Surface -- adiabatic case without tides...*/
			/* Recall that when integrating R, where I should have "v" I still have R components.*/
			y_Riccati_func[0] = 1.0;						//r00	
			y_Riccati_func[1] = -1.0;						//r01
			y_Riccati_func[3] = -(l_func+1.0);				//r11					
			y_Riccati_func[4] = y_Riccati_func[0];			//v0
			y_Riccati_func[5] = y_Riccati_func[1];			//v1		
		}//else (without tides)
		
	}//if(adiabatic_func)
	else
	{
		if(tidesFlag_func)
		{
			y_Riccati_func[0] = 1.0/onePlusZ1;																	//r00
			y_Riccati_func[1] = -(onePlusZ2/onePlusZ1);															//r01
			
			y_Riccati_func[8] = (-4.0*PI*rho_surf_func*surfaceRel_func)/(g_surf_func*onePlusZ1);								//r10
			y_Riccati_func[9] = -1.0 - l_func + (4.0*onePlusZ2*PI*rho_surf_func*surfaceRel_func)/(g_surf_func*onePlusZ1);		//r11
			y_Riccati_func[11] = -(epsilon_t_func*(2.0*l_func + 1.0)*Clmk_func/g_surf_func);					//r13
			
			y_Riccati_func[16] = -(4.0*delAD_V*onePlusZ1 + twoMinus4DelAdV)/(4.0*onePlusZ1);					//r20
			y_Riccati_func[17] = delAD_V + (onePlusZ2*twoMinus4DelAdV)/(4.0*onePlusZ1);							//r21
			y_Riccati_func[18] = 0.25;																			//r22
			
			y_Riccati_func[27] = 1.0;																			//r33
			
		}//if(tidesFlag_func)
		else
		{
			y_Riccati_func[0] = (V_func*(-4.0 + (-1.0 + coeff2)*coeff3 + V_func))/coeff4;			// r00					
			y_Riccati_func[1] = ((1.0 - coeff2*coeff3 + l_func - V_func)*
								 (-4.0 + (-1.0 + coeff2)*coeff3 + V_func) - 
								 4.0*coeff2*(1.0 + coeff2)*wI2*wR2)/coeff4;							// r01
			
			y_Riccati_func[7] = -1.0 - l_func;																			//r11		
			y_Riccati_func[12] = -(delAD_V) + (V_func*(-4.0 + (-1.0 + coeff2)*coeff3 + V_func)*(-0.5 + delAD_V))/coeff4;//r20		
			y_Riccati_func[13] = delAD_V + ((-0.5 + delAD_V)*
											((1.0 - coeff2*coeff3 + l_func - V_func)*
											 (-4.0 + (-1.0 + coeff2)*coeff3 + V_func) - 														
											 4.0*coeff2*(1.0 + coeff2)*wI2*wR2))/coeff4;								//r21
			y_Riccati_func[14] = 0.25;																					//r22			
			y_Riccati_func[3] = (-2.0*(1.0 + coeff2)*V_func*wIwR)/coeff4;												//r03
			y_Riccati_func[4] = (-2.0*(1.0 + l_func + coeff2*(-3.0 - 2.0*coeff3 + l_func) - V_func)*wIwR)/coeff4;		//r04		
			y_Riccati_func[15] = (-2.0*(1.0 + coeff2)*V_func*(-0.5 + delAD_V)*wIwR)/coeff4;								//r23
			y_Riccati_func[16] = (-2.0*(1.0 + l_func + coeff2*(-3.0 - 2.0*coeff3 + l_func) - V_func)*
								  (-0.5 + delAD_V)*wIwR)/coeff4;														//r24			
		}//else without tides
		complete_y_Riccati_usingSimmetryR(y_Riccati_func, rowsR_func);
	}//else non adiabatic

}

void InitialConditionsRij_center(double *y_Riccati_func, int rowsR_func, int adiabatic_func, double l_func, 
								 double c1_func, double omega_r_func, double omega_i_func, int tidesFlag_func)
{	
	int k = 0;
	double wR2 = omega_r_func*omega_r_func, 
	wI2 = omega_i_func*omega_i_func, 
	wIwR = omega_i_func*omega_r_func,
	coeff1 = c1_func * ((wI2 + wR2)*(wI2 + wR2));

	for(k = 0; k < rowsR_func*rowsR_func; k++){y_Riccati_func[k] = 0;}
	if(adiabatic_func)
	{
		if(tidesFlag_func)
		{
			y_Riccati_func[0]  = l_func/(c1_func*wR2);						//r00
			y_Riccati_func[4]  = l_func;									//r11
			y_Riccati_func[8]  = 1.0;										//r22
			y_Riccati_func[9]  = y_Riccati_func[0];							//v0
			y_Riccati_func[10] = y_Riccati_func[1];							//v1
			y_Riccati_func[11] = y_Riccati_func[2];							//v2
		}//if(tidesFlag_func)
		else
		{
			//adiabatic initial conditions
			y_Riccati_func[0] = l_func/(c1_func*wR2);		//r00	
			y_Riccati_func[3] = l_func;						//r11	
			y_Riccati_func[4] = y_Riccati_func[0];			//v0
			y_Riccati_func[5] = y_Riccati_func[1];			//v1					
		}//else without tides
	}//if(adiabatic_func)
	else
	{
		if(tidesFlag_func)
		{
			y_Riccati_func[0] = l_func/(c1_func*wR2);		//r00
			y_Riccati_func[9] = l_func;						//r11
			y_Riccati_func[27] = 1.0;						//r33
		}//if(tidesFlag_func)
		else
		{
			y_Riccati_func[0] = (l_func*(-wI2 + wR2))/coeff1; //r00
			y_Riccati_func[7] = l_func;						  //r11
			y_Riccati_func[3] = (2.0*l_func*wIwR)/coeff1;		//r03
		}//else without tides
		complete_y_Riccati_usingSimmetryR(y_Riccati_func, rowsR_func);
	}//else non adiabatic
	
}


/********************************************************
 ***  At the first time step I am already integrating 
 ***  a permulation of Riccati because of the BCs.
 *** 
 *** - R = 2 x 2 (adiabatic + no tides)
 *** U = {y1, y2}, V = {y3, y4}, to U = {y1, y4}, V = {y2, y3}
 *** 
 *** T = |1 0 0 0|
 ***     |0 0 0 1|
 ***     |0 1 0 0|
 ***     |0 0 1 0|
 *** 
 *** - R = 3 x 3 (adiabatic + tides)
 *** U = {y1, y2, y7}, V = {y3, y4, y8}, to U = {y1, y4, y7}, V = {y2, y3, y8}
 *** 
 *** T = |1 0 0 0 0 0|
 ***     |0 0 0 1 0 0|
 ***     |0 0 1 0 0 0|
 ***     |0 0 0 0 1 0|
 ***     |0 1 0 0 0 0|
 ***     |0 0 0 0 0 1|
 *** 
 *** 
 *** - R = 6 x 6 (non adiabatic + no tides)
 *** U = {y1r, y2r, y5r, y1i, y2i, y5i}, 
 *** V = {y3r, y4r, y6r, y3i, y4i, y6i}, 
 *** to
 *** U = {y1r, y4r, y5r, y1i, y4i, y5i}, 
 *** V = {y2r, y3r, y6r, y2i, y3i, y6i}, 
 *** 
 ***      0 1 2 3 4 5 6 7 8 9 10 11
 *** T = |1 0 0 0 0 0 0 0 0 0 0   0| 0
 ***     |0 0 0 0 0 0 0 1 0 0 0   0| 1
 ***     |0 0 1 0 0 0 0 0 0 0 0   0| 2
 ***     |0 0 0 1 0 0 0 0 0 0 0   0| 3
 ***     |0 0 0 0 0 0 0 0 0 0 1   0| 4
 ***     |0 0 0 0 0 1 0 0 0 0 0   0| 5
 ***     |0 1 0 0 0 0 0 0 0 0 0   0| 6
 ***     |0 0 0 0 0 0 1 0 0 0 0   0| 7
 ***     |0 0 0 0 0 0 0 0 1 0 0   0| 8 
 ***     |0 0 0 0 1 0 0 0 0 0 0   0| 9
 ***     |0 0 0 0 0 0 0 0 0 1 0   0|10
 ***     |0 0 0 0 0 0 0 0 0 0 0   1|11
 *** 
 *** 
//FIXME me. This has been modified.
 *** - R = 8 x 8 (non adiabatic + tides)
 *** U = {y1r, y2r, y5r, y7r, y1i, y2i, y5i, y7i}, 
 *** V = {y3r, y4r, y6r, y8r, y3i, y4i, y6i, y8i}, 
 *** to
 *** U = {y1r, y4r, y5r, y7r, y1i, y4i, y5i, y7i}, 
 *** V = {y2r, y3r, y6r, y8r, y2i, y3i, y6i, y8i}, 
 *** 
 ***      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
 *** T = |1 0 0 0 0 0 0 0 0 0 0   0  0  0  0  0| 0
 ***     |0 0 0 0 0 0 0 0 0 1 0   0  0  0  0  0| 1
 ***     |0 0 1 0 0 0 0 0 0 0 0   0  0  0  0  0| 2
 ***     |0 0 0 1 0 0 0 0 0 0 0   0  0  0  0  0| 3
 ***     |0 0 0 0 1 0 0 0 0 0 0   0  0  0  0  0| 4
 ***     |0 0 0 0 0 0 0 0 0 0 0   0  0  1  0  0| 5
 ***     |0 0 0 0 0 0 1 0 0 0 0   0  0  0  0  0| 6
 ***     |0 0 0 0 0 0 0 1 0 0 0   0  0  0  0  0| 7
 ***     |0 1 0 0 0 0 0 0 0 0 0   0  0  0  0  0| 8
 ***     |0 0 0 0 0 0 0 0 1 0 0   0  0  0  0  0| 9
 ***     |0 0 0 0 0 0 0 0 0 0 1   0  0  0  0  0|10
 ***     |0 0 0 0 0 0 0 0 0 0 0   1  0  0  0  0|11
 ***     |0 0 0 0 0 1 0 0 0 0 0   0  0  0  0  0|12
 ***     |0 0 0 0 0 0 0 0 0 0 0   0  1  0  0  0|13
 ***     |0 0 0 0 0 0 0 0 0 0 0   0  0  0  1  0|14
 ***     |0 0 0 0 0 0 0 0 0 0 0   0  0  0  0  1|15
 *** 
 *** where "r" and "i" denote the real and imaginary part.
 *** 
 *********************************************************/

void pickInitialPermutationMatrix(gsl_matrix *matT_func, gsl_vector* vecPermutInit_BC_func,
								  int rowsT_func, int adiabatic_func, int tidesFlag_func)
{
	int k = 0, j = 0;
	gsl_matrix_set_zero(matT_func);	

	if (adiabatic_func)		
	{
		if (tidesFlag_func)
		{
			for (k = 0; k < rowsT_func; k++)
				for (j = 0; j < rowsT_func; j++)
					if((k == 0 && j == 0)||
					   (k == 1 && j == 3)||
					   (k == 2 && j == 2)||
					   (k == 3 && j == 4)||
					   (k == 4 && j == 1)||
					   (k == 5 && j == 5))
					{
						gsl_vector_set(vecPermutInit_BC_func,k,j);
						gsl_matrix_set(matT_func,k,j, 1.0);
					}				
		}//if (tidesFlag_func) + adiabatic
		else
		{
			for (k = 0; k < rowsT_func; k++)
				for (j = 0; j < rowsT_func; j++)
					if((k == 0 && j == 0)||
					   (k == 1 && j == 3)||
					   (k == 2 && j == 1)||
					   (k == 3 && j == 2))
					{
						gsl_vector_set(vecPermutInit_BC_func,k,j);
						gsl_matrix_set(matT_func,k,j, 1.0);
					}				
		}//else non tides + adiabatic
	}//if (adiabatic_func)
	else
	{
		if(tidesFlag_func)
		{
			for (k=0; k<rowsT_func; k++)
				for (j=0; j<rowsT_func; j++)
					if((k == 0 && j == 1)||(k == 1 && j == 9)||
					   (k == 2 && j == 2)||(k == 3 && j == 3)||
					   (k == 4 && j == 5)||(k == 5 && j == 13)||
					   (k == 6 && j == 6)||(k == 7 && j == 7)||
					   (k == 8 && j == 0)||(k ==  9 && j == 8)||
					   (k == 10 && j == 10)||(k == 11 && j == 11)||
					   (k == 12 && j == 4)||(k == 13 && j == 12)||
					   (k == 14 && j == 14)||(k == 15 && j == 15))
					{
						gsl_vector_set(vecPermutInit_BC_func,k,j);
						gsl_matrix_set(matT_func,k,j, 1.0);
					}							
		}//if(tidesFlag_func) + non adiabatic
		else
		{
			
			for (k=0; k<rowsT_func; k++)
				for (j=0; j<rowsT_func; j++)
					if((k == 0 && j == 0)||(k == 1 && j == 7)||
					   (k == 2 && j == 2)||(k == 3 && j == 3)||
					   (k == 4 && j == 10)||(k == 5 && j == 5)||
					   (k == 6 && j == 1)||(k == 7 && j == 6)||
					   (k == 8 && j == 8)||(k == 9 && j == 4)||
					   (k == 10 && j == 9)||(k == 11 && j == 11))
					{
						gsl_vector_set(vecPermutInit_BC_func,k,j);
						gsl_matrix_set(matT_func,k,j, 1.0);
					}							
		}//else non tides + non adiabatic
	}//else or, non adiabatic
	
}//void pickInitialPermutationMatrix(gsl_matrix *matT_func, int rowsT_func, int adiabatic_func, int tidesFlag_func)



/************************************************************************************
 ****** Initial conditions for the eigenfunctions at the fitting point.
 ****** (see Eq. 20 from Takata & Loffler 2004)
 ****** 
 ****** 
 ****** Adiabatic + no tides
 ****** u0 = y1, u1 = y2, v0 = y3, v1 = y4
 ****** 
 ****** Adiabatic + tides
 ****** u0 = y1, u1 = y2, u2 = y7, v0 = y3, v1 = y4, v2 = y8
 ****** 
 ****** Non Adiabatic + no tides
 ****** u0 = y1, u1 = y2, u2 = y5, v0 = y3, v1 = y4, v2 = y6
 ****** 
 ******  
 ****** Note: For the non adiabatic case, by definition 
 ******(see logbook Chapter 25), R is such that:
 ******  
 ****** Rreal = R00 = R11  
 ****** Rimag = R10 = -R01 
 ******
 ****** where Rij are NxN submatrices of the 2Nx2N R.
 ****** 
 ****** 
 ****** 07/11/2013
 ****** 
 ****** - For the case nonAdiabatic+tides or adiabatic+no tides I am using single value 
 ****** decomposition
 ****** - For the case nonAdiabatic+no tides or adiabatic+tides, I am using the technique below
 ****** 
 ****** 
 ****** 
 ************************************************************************************/

void InitialConditionsUV_fittingPoint(gsl_matrix *matR_diff_func, gsl_vector *vecV_func, int adiabatic_func, int tidesFlag_func)
{
	double v0_init = 0.0, v1_init = 0.0, v2_init = 0.0, v3_init = 0.0, v4_init = 0.0, v5_init = 0.0, v6_init = 0.0, v7_init = 0.0;

	
	if(adiabatic_func)
	{
		if(tidesFlag_func)
		{			
			double v2_init;
			
			v2_init = 100.0;
			v1_init = v2_init*(gsl_matrix_get(matR_diff_func, 1,0)*gsl_matrix_get(matR_diff_func, 0,2) - 
							   gsl_matrix_get(matR_diff_func, 0,0)*gsl_matrix_get(matR_diff_func, 1,2))/
			(gsl_matrix_get(matR_diff_func, 0,0)*gsl_matrix_get(matR_diff_func, 1,1) - 
			 gsl_matrix_get(matR_diff_func, 1,0)*gsl_matrix_get(matR_diff_func, 0,1));

			v0_init = - (gsl_matrix_get(matR_diff_func, 0,1)*v1_init + gsl_matrix_get(matR_diff_func, 0,2)*v2_init)
			/gsl_matrix_get(matR_diff_func, 0,0);
			
			gsl_vector_set(vecV_func,0,v0_init);
			gsl_vector_set(vecV_func,1,v1_init);
			gsl_vector_set(vecV_func,2,v2_init);			
		}//if(tidesFlag_func)
		else
		{
			v1_init = -1.0;
			v0_init = -(gsl_matrix_get(matR_diff_func, 0,1)/gsl_matrix_get(matR_diff_func, 0,0))*v1_init;	
			gsl_vector_set(vecV_func,0,v0_init);
			gsl_vector_set(vecV_func,1,v1_init);			
		}//else adiabatic + non tides	
	}//if(adiabatic_func)
	else
	{		
		if(tidesFlag_func)
		{
			double
			/*R00*/
			r00 = gsl_matrix_get(matR_diff_func, 0, 0),
			r01 = gsl_matrix_get(matR_diff_func, 0, 1),

			r03 = gsl_matrix_get(matR_diff_func, 0, 3),
			
			r10 = gsl_matrix_get(matR_diff_func, 1, 0),
			r11 = gsl_matrix_get(matR_diff_func, 1, 1),
			r12 = gsl_matrix_get(matR_diff_func, 1, 2),
			r13 = gsl_matrix_get(matR_diff_func, 1, 3),
			
			r20 = gsl_matrix_get(matR_diff_func, 2, 0),
			r21 = gsl_matrix_get(matR_diff_func, 2, 1),
			r22 = gsl_matrix_get(matR_diff_func, 2, 2),
			r23 = gsl_matrix_get(matR_diff_func, 2, 3),

			/*R01*/					
			r04 = gsl_matrix_get(matR_diff_func, 0, 4),
			r05 = gsl_matrix_get(matR_diff_func, 0, 5),
			r07 = gsl_matrix_get(matR_diff_func, 0, 7),

			r14 = gsl_matrix_get(matR_diff_func, 1, 4),
			r15 = gsl_matrix_get(matR_diff_func, 1, 5),
			r17 = gsl_matrix_get(matR_diff_func, 1, 7),
			
			r24 = gsl_matrix_get(matR_diff_func, 2, 4),
			r25 = gsl_matrix_get(matR_diff_func, 2, 5),
			r26 = gsl_matrix_get(matR_diff_func, 2, 6),
			r27 = gsl_matrix_get(matR_diff_func, 2, 7),
						
			/*R10 = -R01*/
			r40 = - r04,
			r41 = - r05,
			r43 = - r07,
			
			r50 = - r14,
			r51 = - r15,
			r53 = - r17,

			r60 = - r24,
			r61 = - r25,
			r62 = - r26,
			r63 = - r27,

			/*R11 = R00*/		
			r44 = r00,
			r45 = r01,
			r47 = r03,
			
			r54 = r10,
			r55 = r11,
			r57 = r13,

			r64 = r10,
			r65 = r11,
			r66 = r12,
			r67 = r13,

			V1V0 = r01/r00,
			V3V0 = r03/r00,
			V4V0 = r04/r00,
			V5V0 = r05/r00,
			V7V0 = r07/r00,
			
			V3V1 = (-V3V0*r10 + r13)/(V1V0*r10 - r11),
			V4V1 = (-V4V0*r10 + r14) /(V1V0*r10 - r11),
			V5V1 = (-V5V0*r10 + r15)/(V1V0*r10 - r11),
			V7V1 = (-V7V0*r10 + r17)/(V1V0*r10 - r11),
			
			V3V2 = -((-V3V0*r20 - V1V0*V3V1*r20 + V3V1*r21 + r23)/r22),
			V4V2 = -((-V4V0*r20 - V1V0*V4V1*r20 + V4V1*r21 + r24)/r22),
			V5V2 = -((-V5V0*r20 - V1V0*V5V1*r20 + V5V1*r21 + r25)/r22),
			V6V2 = -(r26/r22),
			V7V2 = -((-V7V0*r20 - V1V0*V7V1*r20 + V7V1*r21 + r27)/r22),
			
			V3V4 = (-V3V0*r40 - V1V0*V3V1*r40 + V3V1*r41 + r43)/(V4V0*r40 + V1V0*V4V1*r40 - V4V1*r41 - r44),
			V5V4 = (-V5V0*r40 - V1V0*V5V1*r40 + V5V1*r41 + r45)/(V4V0*r40 + V1V0*V4V1*r40 - V4V1*r41 - r44),
			V7V4 = (-V7V0*r40 - V1V0*V7V1*r40 + V7V1*r41 + r47)/(V4V0*r40 + V1V0*V4V1*r40 - V4V1*r41 - r44),
			
			
			V3V5 = (-V3V0*r50 - V1V0*V3V1*r50 - V3V4*V4V0*r50 - V1V0*V3V4*V4V1*r50 + 
					V3V1*r51 + V3V4*V4V1*r51 + r53 + V3V4*r54)/
			(V5V0*r50 + V4V0*V5V4*r50 + V1V0*(V5V1 + V4V1*V5V4)*r50 - V5V1*r51 - V4V1*V5V4*r51 - V5V4*r54 - r55),
			
			V7V5 = (-V7V0*r50 - V1V0*V7V1*r50 - V4V0*V7V4*r50 - V1V0*V4V1*V7V4*r50 + 						
					V7V1*r51 + V4V1*V7V4*r51 + V7V4*r54 + r57)/
			(V5V0*r50 + V4V0*V5V4*r50 + V1V0*(V5V1 + V4V1*V5V4)*r50 - V5V1*r51 - V4V1*V5V4*r51 - V5V4*r54 - r55),
			
			
			V3V6 = -(1.0/(V6V2*r62 + r66))*(-V3V0*r60 - V1V0*V3V1*r60 - V3V4*V4V0*r60 - V1V0*V3V4*V4V1*r60 - 												  
											V3V5*V5V0*r60 - V1V0*V3V5*V5V1*r60 - V3V5*V4V0*V5V4*r60 - 												
											V1V0*V3V5*V4V1*V5V4*r60 + V3V1*r61 + V3V4*V4V1*r61 + 
											V3V5*V5V1*r61 + V3V5*V4V1*V5V4*r61 + V3V2*r62 + V3V4*V4V2*r62 + 
											V3V5*V5V2*r62 + V3V5*V4V2*V5V4*r62 + r63 + V3V4*r64 + V3V5*V5V4*r64 + V3V5*r65),
			
			V7V6 = -(1.0/(V6V2*r62 + r66))*(-V7V0*r60 - V1V0*V7V1*r60 - V4V0*V7V4*r60 - V1V0*V4V1*V7V4*r60 - 
											V5V0*V7V5*r60 - V1V0*V5V1*V7V5*r60 - V4V0*V5V4*V7V5*r60 - 
											V1V0*V4V1*V5V4*V7V5*r60 + V7V1*r61 + V4V1*V7V4*r61 + V5V1*V7V5*r61 + 
											V4V1*V5V4*V7V5*r61 + V7V2*r62 + V4V2*V7V4*r62 + V5V2*V7V5*r62 + 
											V4V2*V5V4*V7V5*r62 + V7V4*r64 + V5V4*V7V5*r64 + V7V5*r65 + r67);
			
			
			/*If this part is used, then the IC do not matter, since the tidal eigenfunctions are normalized at the end of the
			 calculation to comply with the boundary conditions*/
			v3_init = 1.0;
			v7_init = 0.0;
			
			v6_init = V3V6*v3_init + V7V6*v7_init;
			v5_init = V3V5*v3_init + V7V5*v7_init;
			v4_init = V3V4*v3_init + V5V4*v5_init + V7V4 * v7_init;
			v2_init = V3V2*v3_init + V4V2*v4_init + V5V2*v5_init + V6V2*v6_init + V7V2*v7_init;
			v1_init = V3V1*v3_init + V4V1*v4_init + V5V1*v5_init + V7V1*v7_init;
			v0_init = -V1V0*v1_init - V3V0*v3_init - V4V0*v4_init - V5V0*v5_init - V7V0*v7_init;

			gsl_vector_set(vecV_func,0,v0_init);
			gsl_vector_set(vecV_func,1,v1_init);
			gsl_vector_set(vecV_func,2,v2_init);
			gsl_vector_set(vecV_func,3,v3_init);
			gsl_vector_set(vecV_func,4,v4_init);
			gsl_vector_set(vecV_func,5,v5_init);
			gsl_vector_set(vecV_func,6,v6_init);
			gsl_vector_set(vecV_func,7,v7_init);
			
		}//if(tidesFlag_func)
		else
		{
			v1_init = -1.0;
			
			double
			/*R00*/
			r00 = gsl_matrix_get(matR_diff_func, 0, 0),
			r01 = gsl_matrix_get(matR_diff_func, 0, 1),
			r02 = gsl_matrix_get(matR_diff_func, 0, 2),
			
			r10 = gsl_matrix_get(matR_diff_func, 1, 0),
			r11 = gsl_matrix_get(matR_diff_func, 1, 1),
			r12 = gsl_matrix_get(matR_diff_func, 1, 2),
			
			r20 = gsl_matrix_get(matR_diff_func, 2, 0),
			r21 = gsl_matrix_get(matR_diff_func, 2, 1),
			r22 = gsl_matrix_get(matR_diff_func, 2, 2),
			
			/*R01*/					
			r03 = gsl_matrix_get(matR_diff_func, 0, 3),
			r04 = gsl_matrix_get(matR_diff_func, 0, 4),
			r05 = gsl_matrix_get(matR_diff_func, 0, 5),
			
			r13 = gsl_matrix_get(matR_diff_func, 1, 3),
			r14 = gsl_matrix_get(matR_diff_func, 1, 4),
			r15 = gsl_matrix_get(matR_diff_func, 1, 5),
			
			r23 = gsl_matrix_get(matR_diff_func, 2, 3),
			r24 = gsl_matrix_get(matR_diff_func, 2, 4),
			r25 = gsl_matrix_get(matR_diff_func, 2, 5),
			
			/*R10 = -R01*/
			r30_diff = - r03,
			r31_diff = - r04,
			r32_diff = - r05,
			
			r40_diff = - r13,
			r41_diff = - r14,
			r42_diff = - r15,
			
			r50_diff = - r23,
			r51_diff = - r24,
			r52_diff = - r25,
			
			/*R11 = R00*/		
			r33_diff = r00,
			r34_diff = r01,
			r35_diff = r02,
			
			r43_diff = r10,
			r44_diff = r11,
			r45_diff = r12,
			
			r53_diff = r20,
			r54_diff = r21,
			r55_diff = r22,
			
			V2V0 = r20/r22,
			V2V1 = r21/r22,
			V2V3 = r23/r22,
			V2V4 = r24/r22,
			V2V5 = r25/r22,
			
			V5det = (-V2V5*r52_diff + r55_diff),
			V5V0 = (r50_diff - V2V0*r52_diff)/V5det,
			V5V1 = (r51_diff - V2V1*r52_diff)/V5det,
			V5V3 = (r53_diff - V2V3*r52_diff)/V5det,
			V5V4 = (r54_diff - V2V4*r52_diff)/V5det,
			
			
			V0det = (r00 - V2V0*r02 + V2V5*V5V0*r02 - V5V0*r05),	
			V0V1 = (r01 - V2V1*r02 + V2V5*V5V1*r02 - V5V1*r05)/V0det,
			V0V3 = (-(V2V3*r02) + V2V5*V5V3*r02 + r03 - V5V3*r05)/V0det,
			V0V4 = (-(V2V4*r02) + V2V5*V5V4*r02 + r04 - V5V4*r05)/V0det,
			
			V3det = (-(V2V3*r32_diff) + V2V5*V5V3*r32_diff + r33_diff - V5V3*r35_diff + V0V3*(-r30_diff + V2V0*r32_diff - V2V5*V5V0*r32_diff + V5V0*r35_diff)),
			V3V1 = (r31_diff - V2V1*r32_diff + V2V5*V5V1*r32_diff - V5V1*r35_diff + V0V1*(-r30_diff + V2V0*r32_diff - V2V5*V5V0*r32_diff + V5V0*r35_diff))/V3det,
			V3V4 = (-(V2V4*r32_diff) + V2V5*V5V4*r32_diff + r34_diff - V5V4*r35_diff + V0V4*(-r30_diff + V2V0*r32_diff - V2V5*V5V0*r32_diff + V5V0*r35_diff))/V3det,
			
			
			V4V1 = (-(V0V1*r40_diff) + V0V3*V3V1*r40_diff + r41_diff + V0V1*V2V0*r42_diff - V2V1*r42_diff - V0V3*V2V0*V3V1*r42_diff + V2V3*V3V1*r42_diff - 
					V0V1*V2V5*V5V0*r42_diff + V0V3*V2V5*V3V1*V5V0*r42_diff + V2V5*V5V1*r42_diff - V2V5*V3V1*V5V3*r42_diff - V3V1*r43_diff + 
					V0V1*V5V0*r45_diff - V0V3*V3V1*V5V0*r45_diff - V5V1*r45_diff + V3V1*V5V3*r45_diff)/
			(-(V0V4*r40_diff) + V0V3*V3V4*r40_diff + V0V4*V2V0*r42_diff - V2V4*r42_diff - V0V3*V2V0*V3V4*r42_diff + V2V3*V3V4*r42_diff - 
			 V0V4*V2V5*V5V0*r42_diff + V0V3*V2V5*V3V4*V5V0*r42_diff - V2V5*V3V4*V5V3*r42_diff + V2V5*V5V4*r42_diff - V3V4*r43_diff + 
			 r44_diff + V0V4*V5V0*r45_diff - V0V3*V3V4*V5V0*r45_diff + V3V4*V5V3*r45_diff - V5V4*r45_diff);
			
			v4_init = -V4V1*v1_init;
			v3_init = -V3V1*v1_init - V3V4*v4_init;
			v0_init = -V0V1*v1_init - V0V3*v3_init - V0V4*v4_init;
			v5_init = -V5V0*v0_init - V5V1*v1_init - V5V3*v3_init - V5V4*v4_init;
			v2_init = -V2V0*v0_init - V2V1*v1_init - V2V3*v3_init - V2V4*v4_init - V2V5*v5_init;
			
			
			gsl_vector_set(vecV_func,0,v0_init);
			gsl_vector_set(vecV_func,1,v1_init);
			gsl_vector_set(vecV_func,2,v2_init);
			gsl_vector_set(vecV_func,3,v3_init);
			gsl_vector_set(vecV_func,4,v4_init);
			gsl_vector_set(vecV_func,5,v5_init);
		}//else no tides
			
	}//else non adiabatic



}//void InitialConditionsUV_fittingPoint