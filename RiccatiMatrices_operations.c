/*************************************************************************************
 ************ 
 ************	Riccati matrices configuration:
 ************ 
 ************ adiabatic, no tides:
 ************ 
 ************ U = y1, y2
 ************ V = y3, y4
 ************ 
 ************ adiabatic, tides:
 ************ 
 ************ U = y1, y2, y7
 ************ V = y3, y4, y8
 ************ 
 ************ non adiabatic, no tides:
 ************ 
 ************ U = y1r, y2r, y5r, y1i, y2i,y5i
 ************ V = y3r, y4r, y6r, y3i, y4i, y6i 
 ************ 
 ************ non adiabatic, tides:
 ************ 
 ************ U = y1r, y2r, y5r, y7r, y1i, y2i, y5i, y7i
 ************ V = y3r, y4r, y6r, y8r, y3i, y4i, y6i, y8i 
 ************ 
 **************************************************************************************/
#include <iomanip>
#include <iostream>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include "RiccatiMatrices_operations.h"
#include "matrix_operations.h"
#include "IOfiles.h"

using namespace std;

/*************************************************************************************
 ************ 
 ************	The routines described below are used to perform operations on the 
 ************	riccati matrices. The matrices involved are:
 ************ 
 ************ - R = of size (rowsR, rowsR)
 ************ - M = of size (2*rowsR, 2*rowsR) which contains the coefficients of the differential
 ************  equations for y_i
 ************ - A, B, C, D, which are submatrices of M of size (rowsR, rowsR)
 ************ 
 **************************************************************************************/


/*************************************************************************************
 ************ 
 ************	This routine creates the Riccati matrices A, B, C, D.
 ************	In input are the various stellar quantities. Note that the stellar
 ************	pulsation equations have been modified such that the integration 
 ************   variables are y_i'= y_i/r^(l-2)
 ************ 
 **************************************************************************************/

void create_Riccati_ABCD_new(double csi_func, double Vg_func,double Astar_func, double U_func,double l_func,
							 double omega_r_func, double omega_i_func, double Vt_func, double V_func, double del_func, double delad_func, 							
							 double ks_func, double c1_func, double c2_func, double c3_func, double c4_func, double epsAd_func, double dlnLR_dlnr_func, 
							 double epsS_func, gsl_matrix *A_func, gsl_matrix *B_func, gsl_matrix *C_func, gsl_matrix *D_func, int rowsRic_func, 
							 int inWard_outWard_func, int nonAdiabatic_func, int tidesFlag_func, int WD_tides_flag_func)
{
	
	double wR2 = omega_r_func*omega_r_func,
	wI2 = omega_i_func*omega_i_func, 
	wRwI = omega_i_func*omega_r_func, 
	llp1 = l_func*(1.0 + l_func), 
	c2V = c2_func*V_func,
	delMdelAd = del_func - delad_func,
	c1wR2 = c1_func*wR2,
	deladV = delad_func*V_func,
	delV = del_func*V_func,
	wI2pwR2_2 = (wI2 + wR2)*(wI2 + wR2);
	int i = 0, j = 0;
    double changeVariables = csi_func;

    //FIXME, change csti with P everywhere!!!!!
    if(WD_tides_flag_func){changeVariables = V_func;}
    
    
	gsl_matrix_set_zero(A_func);
	gsl_matrix_set_zero(B_func);
	gsl_matrix_set_zero(C_func);
	gsl_matrix_set_zero(D_func);

	
    
	/*In the purely adiabatic case omega_i = 0, and the 2x2 matrices below
	 boil down to the purely adiabatic case...*/
	gsl_matrix_set(A_func,0,0,(-1.0 - l_func + Vg_func)/changeVariables);
	gsl_matrix_set(A_func,0,1,(-Vg_func + (llp1*(-wI2 + wR2)/(c1_func*wI2pwR2_2)))/changeVariables);		
	gsl_matrix_set(A_func,1,0,(-Astar_func + c1_func*(-wI2 + wR2))/changeVariables);
	gsl_matrix_set(A_func,1,1,(3.0 + Astar_func - l_func - U_func)/changeVariables);
	
	gsl_matrix_set(B_func,0,0,Vg_func/changeVariables);
	gsl_matrix_set(B_func,0,1,0.0);
	gsl_matrix_set(B_func,1,0,-Astar_func/changeVariables);
	gsl_matrix_set(B_func,1,1,0.0);		

	
	gsl_matrix_set(C_func,0,0,0.0);
	gsl_matrix_set(C_func,0,1,0.0);
	gsl_matrix_set(C_func,1,0,(Astar_func*U_func)/changeVariables);
	gsl_matrix_set(C_func,1,1,(U_func*Vg_func)/changeVariables);		
	
	gsl_matrix_set(D_func,0,0,(3.0 - l_func - U_func)/changeVariables);
	gsl_matrix_set(D_func,0,1, 1.0/changeVariables);
	gsl_matrix_set(D_func,1,0,(llp1 - U_func*Vg_func)/changeVariables);
	gsl_matrix_set(D_func,1,1,(2.0 - l_func - U_func)/changeVariables);			

	if(nonAdiabatic_func == 0 && tidesFlag_func)
	{
		gsl_matrix_set(A_func,2,2,(2.0 - l_func)/changeVariables);
		gsl_matrix_set(D_func,2,2,(2.0 - l_func)/changeVariables);		
	}
	
	if(nonAdiabatic_func)
	{
		if(tidesFlag_func)
		{
			gsl_matrix_set(A_func,0,2,Vt_func/changeVariables);
			gsl_matrix_set(A_func,1,2,Vt_func/changeVariables);
			gsl_matrix_set(A_func,2,0,V_func*(c2_func + 4.0*del_func + delad_func*(-4.0 + U_func - c1wR2))/changeVariables);			
			gsl_matrix_set(A_func,2,1,(-c2V - (delMdelAd*llp1*V_func)/(c1wR2))/changeVariables);			
			gsl_matrix_set(A_func,2,2,(2.0 - l_func - del_func*(-4.0 + ks_func)*V_func)/changeVariables);
			gsl_matrix_set(A_func,3,3,(2.0 - l_func)/changeVariables);
			
			gsl_matrix_set(B_func,2, 0, c2V/changeVariables);			
			gsl_matrix_set(B_func,2, 1, deladV/changeVariables);			
			gsl_matrix_set(B_func,2, 2, -delV/changeVariables);
						
			gsl_matrix_set(C_func,1,2,-U_func*Vt_func/changeVariables);
			gsl_matrix_set(C_func,2,0,((delad_func*llp1 - del_func*(llp1 + c3_func*epsAd_func*V_func))/del_func)/changeVariables);
			gsl_matrix_set(C_func,2,1,((-(c1wR2*delad_func*llp1) + c3_func*del_func*(llp1 + c1wR2*epsAd_func*V_func))/(c1wR2*del_func))/changeVariables);
			gsl_matrix_set(C_func,2,2,(-((llp1 - c3_func*delV*epsS_func)/delV))/changeVariables);
			gsl_matrix_set(C_func,2,6,c4_func*omega_r_func/changeVariables);
			
			gsl_matrix_set(D_func,2,0,((delad_func*llp1 - c3_func*del_func*V_func*epsAd_func)/del_func)/changeVariables);
			gsl_matrix_set(D_func,2,2,(2 - dlnLR_dlnr_func - l_func)/changeVariables);
			gsl_matrix_set(D_func,3,3,(2.0 - l_func)/changeVariables);
		}//if(tidesFlag_func)
		else
		{
			gsl_matrix_set(A_func,0,2,Vt_func/changeVariables);
			gsl_matrix_set(A_func,0,4,((2.0*llp1*wRwI)/(c1_func*wI2pwR2_2))/changeVariables);			
			gsl_matrix_set(A_func,1,2,(Vt_func)/changeVariables);
			gsl_matrix_set(A_func,1,3,(-2.0*c1_func*wRwI)/changeVariables);
			gsl_matrix_set(A_func,2,0,(V_func*(c2_func + 4.0*del_func + delad_func*(-4.0 + U_func + c1_func*(wI2 - wR2))))/changeVariables);
			gsl_matrix_set(A_func,2,1,(-c2V + (delMdelAd*l_func*(1.0 + l_func)*V_func*(wI2 - wR2))/
									   (c1_func*wI2pwR2_2))/changeVariables);
			gsl_matrix_set(A_func,2,2,(2.0 - l_func + del_func*(4.0 - ks_func)*V_func)/changeVariables);
			gsl_matrix_set(A_func,2,3,(2.0*c1_func*deladV*wRwI)/changeVariables);
			gsl_matrix_set(A_func,2,4,((-2.0*delMdelAd*llp1*V_func*wRwI)/(c1_func*wI2pwR2_2))/changeVariables);
			
			gsl_matrix_set(B_func,2,0,c2V/changeVariables);
			gsl_matrix_set(B_func,2,1,deladV/changeVariables);
			gsl_matrix_set(B_func,2,2,-delV/changeVariables);
			
			
			gsl_matrix_set(C_func,1,2,(-(U_func*Vt_func))/changeVariables);	
			gsl_matrix_set(C_func,2,0,(-((delMdelAd*llp1)/del_func) - c3_func*epsAd_func*V_func)/changeVariables);
			gsl_matrix_set(C_func,2,1,((-(delad_func*llp1) + c3_func*del_func*epsAd_func*V_func)/del_func - 
									   (c3_func*llp1*(wI2 - wR2))/(c1_func*wI2pwR2_2))/changeVariables);
			
			gsl_matrix_set(C_func,2,2,((-(llp1) + c3_func*del_func*epsS_func*V_func + 
										c4_func*delV*omega_i_func)/delV)/changeVariables);
			
			gsl_matrix_set(C_func,2,3,0.0);
			gsl_matrix_set(C_func,2,4,((2.0*c3_func*llp1*wRwI)/(c1_func*wI2pwR2_2))/changeVariables);		
			gsl_matrix_set(C_func,2,5,(c4_func*omega_r_func)/changeVariables);	
			
			
			gsl_matrix_set(D_func,2,0,((delad_func*llp1)/del_func - c3_func*epsAd_func*V_func)/changeVariables);
			gsl_matrix_set(D_func,2,2,(2.0 - dlnLR_dlnr_func - l_func)/changeVariables);
		}
		
//		/*Finish A, B, C, D using simmetry properties*/
		for(i = rowsRic_func/2; i<rowsRic_func; i++)
			for(j = rowsRic_func/2; j<rowsRic_func; j++)
			{
				gsl_matrix_set(A_func,i,j, gsl_matrix_get(A_func, i-rowsRic_func/2, j-rowsRic_func/2));	
				gsl_matrix_set(B_func,i,j, gsl_matrix_get(B_func, i-rowsRic_func/2, j-rowsRic_func/2));	
				gsl_matrix_set(C_func,i,j, gsl_matrix_get(C_func, i-rowsRic_func/2, j-rowsRic_func/2));	
				gsl_matrix_set(D_func,i,j, gsl_matrix_get(D_func, i-rowsRic_func/2, j-rowsRic_func/2));	
				
				gsl_matrix_set(A_func,i,j-rowsRic_func/2, -gsl_matrix_get(A_func, i-rowsRic_func/2, j));	
				gsl_matrix_set(B_func,i,j-rowsRic_func/2, -gsl_matrix_get(B_func, i-rowsRic_func/2, j));	
				gsl_matrix_set(C_func,i,j-rowsRic_func/2, -gsl_matrix_get(C_func, i-rowsRic_func/2, j));	
				gsl_matrix_set(D_func,i,j-rowsRic_func/2, -gsl_matrix_get(D_func, i-rowsRic_func/2, j));					
			}		
	}//	if(nonAdiabatic_func)
	
}//create_Riccati_ABCD_new


void permute_Riccati_ABCD(double csi_func, double Vg_func,double Astar_func, double U_func,double c1_func, double l_func,
						  double omega_r_func, double omega_i_func, double Vt_func, double V_func, double del_func, double delad_func, 							
						  double ks_func, double c2_func, double c3_func, double c4_func, double epsAd_func, double dlnLR_dlnr_func, 
						  double epsS_func, 
						  gsl_matrix *A_func, gsl_matrix *B_func, gsl_matrix *C_func, gsl_matrix *D_func, gsl_matrix *M_func, 
						  gsl_matrix *T_func,gsl_matrix *matDummySizeT_func, gsl_matrix *matDummySizeT_2_func, 
						  gsl_vector *vecDummySizeT2_func, gsl_permutation *permDummySizeT_func, int rowsRic_func, 
						  int inWard_outWard_func, int nonAdiabatic_func, int tidesFlag_func, int WD_tides_flag_func)
{
	
	int rowsT_func = 2*rowsRic_func, i = 0, j = 0;
	
	/*** this is the form of the original riccati matrix 
	 *** 
	 *** - R = 2 x 2
	 *** U = {y1, y2}, V = {y3, y4}
	 *** 
	 *** - R = 6 x 6
	 *** U = {y1r, y2r, y5r, y1i, y2i, y5i}, 
	 *** V = {y3r, y4r, y6r, y3i, y4i, y6i}, 
	 *** 
	 *** With y divided by r^(l-2)	
	 ***/
	
		create_Riccati_ABCD_new(csi_func, Vg_func, Astar_func, U_func, l_func,
								omega_r_func, omega_i_func, Vt_func, V_func, del_func, delad_func, 							
								ks_func, c1_func, c2_func, c3_func, c4_func, epsAd_func, dlnLR_dlnr_func, 
								epsS_func, A_func, B_func, C_func, D_func, rowsRic_func, inWard_outWard_func,
								nonAdiabatic_func, tidesFlag_func, WD_tides_flag_func);

	/*storing the submatrices A, B, C, and D in M */
	for (i=0; i<rowsRic_func; i++)
		for (j=0; j<rowsRic_func; j++)
		{
			gsl_matrix_set(M_func, i, j, gsl_matrix_get(A_func, i, j));
			gsl_matrix_set(M_func, i, j+rowsRic_func, gsl_matrix_get(B_func, i, j));
			gsl_matrix_set(M_func, i+rowsRic_func, j, gsl_matrix_get(C_func, i, j));
			gsl_matrix_set(M_func, i+rowsRic_func, j+rowsRic_func, gsl_matrix_get(D_func, i, j));			
		}	
	/********************************************************
	 ***** 
	 *****    Applying Eq (27) from Takata Loffler 2004
	 *****		     2004PASJ...56..645T
	 ***** 
	 *****			   M' = TM(T^-1)
	 *********************************************************/
	
	/*invert matrix T and write it on a dummy matrix*/
	InverseMatrix(T_func, matDummySizeT_func, rowsT_func, permDummySizeT_func, matDummySizeT_2_func);
	
	/* Overwrite MT^-1 on the dummy matrix*/
	matrix_mul_rows_by_cols(M_func, matDummySizeT_func, matDummySizeT_func, rowsT_func, vecDummySizeT2_func);
	
	/* Overwrite TMT^-1 on the dummy matrix*/
	matrix_mul_rows_by_cols(T_func, matDummySizeT_func, matDummySizeT_func, rowsT_func, vecDummySizeT2_func);
	
	/*Now the dummy matrix == M', and I can exctract my new A', B', C', D' */
	extract_submatrices(matDummySizeT_func, A_func, B_func, C_func, D_func, rowsRic_func);	
}//permute_Riccati_ABCD

void permute_Riccati_matrix_R(gsl_matrix *T_func, gsl_matrix *R_func, 
							  gsl_matrix *T00_func,gsl_matrix *T01_func,gsl_matrix *T10_func,gsl_matrix *T11_func, 
							  gsl_matrix *matRperm_func, gsl_matrix *matDummySizeR_func, gsl_matrix *matDummySizeR_2_func,
							  gsl_matrix *matDummySizeR_3_func, gsl_vector *vecDummySizeR2_func, gsl_permutation *permDummySizeR_func, 
							  int rowsRic_func)
{
	/********************************************************
	 ***** 
	 *****    Applying Eq (25) from Takata Loffler 2004
	 *****		     2004PASJ...56..645T
	 ***** 
	 *****	  R' = (T00*R + T01)(T10*R + T11)^-1	 
	 *********************************************************/
	double determinant = 0.0;

	extract_submatrices(T_func, T00_func, T01_func, T10_func, T11_func, rowsRic_func);
		
	/* Overwrite T10*R on matDummySizeR_func*/
	matrix_mul_rows_by_cols(T10_func, R_func, matDummySizeR_func, rowsRic_func, vecDummySizeR2_func);
	/* Now, T10*R = matDummySizeR_func*/

	/* Overwrite T10*R + T11 on matDummySizeR_func*/
	matrix_sum(matDummySizeR_func, T11_func, matDummySizeR_func, rowsRic_func);
	/* Now, T10*R + T11 = matDummySizeR_func*/
	
	
	determinant = determinantNbyN(matDummySizeR_func, rowsRic_func, permDummySizeR_func, matDummySizeR_2_func);

	/* check determinant */
	if (determinant == 0.0)
	{
		errorMessageAndExit("RiccatiMatrices_operations.c (permute_Riccati_matrix_R)", 
							"SINGULAR MATRIX...try increasing Riccati limit");
	}
	else
	{
		/* Overwrite (T10*R + T11)^-1 on matDummySizeR_func*/
		InverseMatrix(matDummySizeR_func, matDummySizeR_3_func, rowsRic_func, permDummySizeR_func, matDummySizeR_2_func);
		/* Now, (T10*R + T11)^-1 = matDummySizeR_func*/
		
		/* Overwrite T00*R on matRperm_func*/
		matrix_mul_rows_by_cols(T00_func, R_func, matRperm_func, rowsRic_func, vecDummySizeR2_func);
		/* Now, T00*R = matRperm_func*/

		/* Overwrite T00*R + T01 on matRperm_func*/
		matrix_sum(matRperm_func, T01_func, matRperm_func, rowsRic_func);	
		/* Now, T00*R + T01 = matRperm_func*/		

		/* Overwrite(T00*R + T01)(T10*R + T11)^-1 on matRperm_func */	
		matrix_mul_rows_by_cols(matRperm_func, matDummySizeR_3_func, matRperm_func, rowsRic_func, vecDummySizeR2_func);
	}
}//permute_Riccati_matrix_R

void create_Riccati_eqs(gsl_matrix *A_func, gsl_matrix *B_func, gsl_matrix *C_func, gsl_matrix *D_func, 
						gsl_matrix *R_func, gsl_matrix *dR_func, int rowsRic_func, gsl_matrix *matDummySizeR_func, 
						gsl_matrix *matDummySizeR_2_func, gsl_vector *vecDummySizeR2_func)
{
	
	/********************************************************
	 ***** 
	 *****       Using Eq. (13) from Takata Loffler 2004
	 *****				2004PASJ...56..645T
	 *****           dR/dx = B + AR - RD - RCR
	 ***** 
	 *********************************************************/
	/* AR and overwrite matDummySizeR_func --> AR = matDummySizeR_func */
	matrix_mul_rows_by_cols(A_func, R_func, matDummySizeR_func, rowsRic_func, vecDummySizeR2_func);

	/* B + AR and overwrite matDummySizeR_func --> B + AR = matDummySizeR_func*/
	matrix_sum (B_func, matDummySizeR_func, matDummySizeR_func, rowsRic_func);

	/* RD and overwrite matDummySizeR_2_func --> RD = matDummySizeR_2_func*/
	matrix_mul_rows_by_cols(R_func, D_func, matDummySizeR_2_func, rowsRic_func, vecDummySizeR2_func);
	
	/* B + AR - RD === (matDummySizeR_func - matDummySizeR_2_func) and overwrite matDummySizeR_func
	 --> B + AR - RD = matDummySizeR_func */	
	matrix_diff (matDummySizeR_func, matDummySizeR_2_func, matDummySizeR_func, rowsRic_func);	
	
	/* RC and overwrite matDummySizeR_2_func --> RC = matDummySizeR_2_func*/
	matrix_mul_rows_by_cols(R_func, C_func, matDummySizeR_2_func, rowsRic_func, vecDummySizeR2_func);

	/* RCR = (matDummySizeR_2_func*R) and overwrite matDummySizeR_2_func --> RCR = matDummySizeR_2_func*/
	matrix_mul_rows_by_cols(matDummySizeR_2_func, R_func, matDummySizeR_2_func, rowsRic_func, vecDummySizeR2_func);

	/* B + AR - RD - RCR = (matDummySizeR_func - matDummySizeR_2_func) and write on dR_func */
	matrix_diff (matDummySizeR_func, matDummySizeR_2_func, dR_func, rowsRic_func);
}//create_Riccati_eqs

double calc_detRiccati_condition(gsl_matrix *T_func, gsl_matrix *T00_func, gsl_matrix *T01_func, 
								 gsl_matrix *T10_func,gsl_matrix *T11_func, gsl_matrix *R_func, 
								 gsl_matrix *matDummySizeR_func, gsl_matrix *matDummySizeR_2_func, 
								 gsl_vector *vecDummySizeR2_func, gsl_permutation *permDummySizeR_func, 
								 int rowsRic_func)
{
	double determinant = 0.0;
	/********************************************************
	 ***** 
	 *****    performing check Eq (26) from Takata Loffler 2004
	 *****		     2004PASJ...56..645T
	 ***** 
	 *****	           det(T10*R + T11) != 0
	 *********************************************************/
	extract_submatrices(T_func, T00_func, T01_func, T10_func, T11_func, rowsRic_func);

	matrix_mul_rows_by_cols(T10_func, R_func, matDummySizeR_func, rowsRic_func, vecDummySizeR2_func);											
	
	matrix_sum(matDummySizeR_func, T11_func, matDummySizeR_func, rowsRic_func);											

	determinant = determinantNbyN(matDummySizeR_func, rowsRic_func, permDummySizeR_func, matDummySizeR_2_func);

	return determinant;
}//calc_detRiccati_condition

void back_to_original_Riccati_R(gsl_matrix *T_func, gsl_matrix *T00_func, gsl_matrix *T01_func, 
								gsl_matrix *T10_func,gsl_matrix *T11_func, gsl_matrix *R_func, 
								gsl_matrix *R_original_func, gsl_matrix *matDummySizeR_func, 
								gsl_matrix *matDummySizeR_2_func, gsl_matrix *matDummySizeR_3_func,
								gsl_vector *vecDummySizeR2_func, gsl_permutation *permDummySizeR_func, 
								int rowsRic_func)
{
	
	/********************************************************
	 ***** 
	 *****    Inverting Eq (25) from Takata Loffler 2004
	 *****		     2004PASJ...56..645T
	 ***** 
	 *****	     R = (R'T10 - T00)^-1(T01 - R'T11)
	 *********************************************************/
	double determinant = 0.0;
	extract_submatrices(T_func, T00_func, T01_func, T10_func, T11_func, rowsRic_func);
	
	/* Overwrite R'T11 on matDummySizeR_func*/
	matrix_mul_rows_by_cols(R_func, T11_func, matDummySizeR_func, rowsRic_func, vecDummySizeR2_func);											
			
	/* Overwrite T01 - R'T11 on matDummySizeR_func*/
	matrix_diff(T01_func, matDummySizeR_func, matDummySizeR_func, rowsRic_func);											
	
	/* Overwrite R'T10 on R_original_func*/
	matrix_mul_rows_by_cols(R_func, T10_func, R_original_func, rowsRic_func, vecDummySizeR2_func);											
	
	/* Overwrite R'T10 - T00 on R_original_func*/
	matrix_diff(R_original_func, T00_func, R_original_func, rowsRic_func);	
	
	/*check on the determinant*/
	determinant = determinantNbyN(R_original_func, rowsRic_func, permDummySizeR_func, matDummySizeR_2_func);

	if (determinant == 0.0)
	{
		errorMessageAndExit("RiccatiMatrices_operations.c (back_to_original_Riccati_R)", 
							"SINGULAR MATRIX...try increasing Riccati limit or decrease accuracy requirement");
	}
	else
	{
		/* Overwrite (R'T10 - T00)^-1 on matDummySizeR_3_func*/
		InverseMatrix(R_original_func, matDummySizeR_3_func, rowsRic_func, permDummySizeR_func, matDummySizeR_2_func);
		matrix_mul_rows_by_cols(matDummySizeR_3_func, matDummySizeR_func, R_original_func, rowsRic_func, vecDummySizeR2_func);											
	}
}//back_to_original_Riccati_R


/****************************************************************************************
 ***** 
 *****		FINDING THE PERMUTATION WITH THE MINIMUM ||R||
 ***** 
 ***** - adiabatic (with or without tides): I look for all possible permutations.
 ***** 
 ***** - non adiabatic: the number of permutations to perform sould be T = 12! For R = 6x6.
 ***** To avoid an infinite running time, I perform all the possible permutations on the 
 ***** REAL part of R (3x3 ==> # permutations 6! = 720)
 ***** 
 ***** - non adiabatic + tides: I perform all the possible permutations on the 
 ***** REAL part of R and I keep y7 and y8 fixed.
 ***** 	 
 ******************************************************************/
void find_permutation_with_minNorm(gsl_matrix *matT_func, gsl_matrix *matT00_func, gsl_matrix *matT01_func, 
								   gsl_matrix *matT10_func, gsl_matrix *matT11_func, gsl_matrix *matR_func, 
								   gsl_matrix *matDummySizeR_func, gsl_matrix *matDummySizeR_2_func, 
								   gsl_matrix *matDummySizeR_3_func, gsl_vector *vecPermutIndices_func, 
								   gsl_vector *vecDummySizeR2_func, gsl_vector *vecDummySizeT_func, 
								   gsl_permutation *permDummySizeR_func, int rowsRic_func, int &PermutationMinNormNumber_func, 
								   int adiabatic_func, int tidesFlag_func)
{
	double det_T10R_T11_s = 0.0, normR = 0.0, normR_min = 0.0;
	int normR_counts = 0, k = 0, i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0, i6 = 0, m = 0, n = 0, 
	rowsT_func = 2 * rowsRic_func, allowedPermutations = 0;

	if(adiabatic_func)
	{
		if(tidesFlag_func)
			allowedPermutations = 720;  //for a 6 x 6 size T (6 factorial);		
		else
			allowedPermutations = 24;  //for a 4 x 4 size T (4 factorial);		
	}//else adiabatic 
	else
	{
		if(tidesFlag_func)
			allowedPermutations = 5040;
		else
			allowedPermutations = 720;
	}
	gsl_matrix *mat_normR_func =  gsl_matrix_alloc(allowedPermutations, rowsT_func+1+1);
	gsl_matrix *matRperm_func  =  gsl_matrix_alloc(rowsRic_func, rowsRic_func);
	
	gsl_matrix_set_zero(mat_normR_func);
	
	if(adiabatic_func)		
	{
		if(tidesFlag_func)
		{
			for (i1 = 0; i1 < rowsT_func ; i1++)
			for (i2 = 0; i2 < rowsT_func; i2++)
			if (i2 != i1)		
			for(i3 = 0; i3 < rowsT_func; i3++)
			if (i3 != i1 && i3 != i2)		
			for(i4 = 0; i4 < rowsT_func; i4++)
			if (i4 != i1 && i4 != i2 && i4 != i3)		
			for(i5 = 0; i5 < rowsT_func; i5++)
			if (i5 != i1 && i5 != i2 && i5 != i3 && i5 != i4)		
			for(i6 = 0; i6 < rowsT_func; i6++)
			if (i6 != i1 && i6 != i2 && i6 != i3 && i6 != i4 && i6 != i5)
			{
				gsl_vector_set(vecDummySizeT_func,0,i1);
				gsl_vector_set(vecDummySizeT_func,1,i2);
				gsl_vector_set(vecDummySizeT_func,2,i3);
				gsl_vector_set(vecDummySizeT_func,3,i4);
				gsl_vector_set(vecDummySizeT_func,4,i5);
				gsl_vector_set(vecDummySizeT_func,5,i6);

				gsl_matrix_set_zero(matT_func);
				
				for (m = 0; m<rowsT_func; m++)
					gsl_matrix_set(matT_func,m, static_cast<int>(gsl_vector_get(vecDummySizeT_func,m)), 1.0);
				
				/********************************************************
				 ***** 
				 *****    performing check Eq (26) from Takata Loffler 2004
				 *****		     2004PASJ...56..645T
				 ***** 
				 *****	           det(T10*R + T11) != 0
				 *********************************************************/
				det_T10R_T11_s = calc_detRiccati_condition(matT_func, matT00_func, matT01_func, matT10_func,
														   matT11_func, matR_func, matDummySizeR_func, 
														   matDummySizeR_2_func, vecDummySizeR2_func, 
														   permDummySizeR_func, rowsRic_func);
				
				/*if the condition on the determinant is satisfied, then permute riccati*/
				if(det_T10R_T11_s != 0.0)
				{							
					permute_Riccati_matrix_R(matT_func, matR_func, matT00_func, matT01_func, matT10_func, matT11_func, 
											 matRperm_func, matDummySizeR_func, matDummySizeR_2_func, matDummySizeR_3_func, 
											 vecDummySizeR2_func, permDummySizeR_func, rowsRic_func);
					
					/* and now I have a permutation of R in matRperm_func. Calculate the norm of it, and 
					 store it together with the indices of the permutation. THis will be used later to find
					 the permutation with the minimum norm.*/						
					
					normR = 0.0;
					for (m=0; m<rowsRic_func; m++)
						for(n=0; n<rowsRic_func; n++){normR = normR + pow(gsl_matrix_get(matRperm_func,m,n), 2.0);}
					normR = sqrt(normR);
					
					gsl_matrix_set(mat_normR_func, k, 0,normR);
					
					/*store all permutations and corresponding ||R||*/
					for (m = 1; m<1+rowsT_func; m++)
						gsl_matrix_set(mat_normR_func,k, m, gsl_vector_get(vecDummySizeT_func,m-1));
					
					k++;
				}//if(det_T10R_T11_s != 0.0)											
			}//if (i6 != i1 && i6 != i2 && i6 != i3 && i6 != i4 && i6 != i5)
		}//if(tidesFlag_func)
		else
		{
			/* try each permutation */
			for (i1 = 0; i1 <rowsT_func; i1++)
				for (i2 = 0; i2 <rowsT_func; i2++)
					if (i2 != i1)		
						for (i3 = 0; i3 <rowsT_func; i3++)
							if (i3 != i1 && i3 != i2) 
								for (i4 = 0; i4 <rowsT_func; i4++)
									if (i4 != i1 && i4 != i2 && i4 != i3)									
									{
										gsl_vector_set(vecDummySizeT_func,0,i1);
										gsl_vector_set(vecDummySizeT_func,1,i2);
										gsl_vector_set(vecDummySizeT_func,2,i3);
										gsl_vector_set(vecDummySizeT_func,3,i4);
										gsl_matrix_set_zero(matT_func);
										
										for (m = 0; m<rowsT_func; m++)
											gsl_matrix_set(matT_func,m, static_cast<int>(gsl_vector_get(vecDummySizeT_func,m)), 1.0);
										
										/********************************************************
										 ***** 
										 *****    performing check Eq (26) from Takata Loffler 2004
										 *****		     2004PASJ...56..645T
										 ***** 
										 *****	           det(T10*R + T11) != 0
										 *********************************************************/
										det_T10R_T11_s = calc_detRiccati_condition(matT_func, matT00_func, matT01_func, matT10_func,
																				   matT11_func, matR_func, matDummySizeR_func, 
																				   matDummySizeR_2_func, vecDummySizeR2_func, 
																				   permDummySizeR_func, rowsRic_func);
										
										/*if the condition on the determinant is satisfied, then permute riccati*/
										if(det_T10R_T11_s != 0.0)
										{							
											permute_Riccati_matrix_R(matT_func, matR_func, matT00_func, matT01_func, matT10_func, matT11_func, 
																	 matRperm_func, matDummySizeR_func, matDummySizeR_2_func, matDummySizeR_3_func, 
																	 vecDummySizeR2_func, permDummySizeR_func, rowsRic_func);
											
											/* and now I have a permutation of R in matRperm_func. Calculate the norm of it, and 
											 store it together with the indices of the permutation. THis will be used later to find
											 the permutation with the minimum norm.*/						
											
											normR = 0.0;
											for (m=0; m<rowsRic_func; m++)
												for(n=0; n<rowsRic_func; n++){normR = normR + pow(gsl_matrix_get(matRperm_func,m,n), 2.0);}
											normR = sqrt(normR);
											
											gsl_matrix_set(mat_normR_func, k, 0,normR);
											
											/*store all permutations and corresponding ||R||*/
											for (m = 1; m<1+rowsT_func; m++)
												gsl_matrix_set(mat_normR_func,k, m, gsl_vector_get(vecDummySizeT_func,m-1));
											
											k++;
										}//if(det_T10R_T11_s != 0.0)											
									}//if (i4 != i1 && i4 != i2 && i4 != i3)					
		}//else no tides + adiabatic
	}//	if(adiabatic_func)		
	else if (adiabatic_func == 0)
	{
		int permutationCounter = 0;
		int i7 = 0, i8 = 0, i9 = 0, 		
		index1 = 0, index2 = 0, index3 = 0, index4 = 0, index5 = 0, index6 = 0, 
		index7 = 0, index8 = 0, index9 = 0, index10 = 0, index11 = 0, index12 = 0; 
		
		
		if(tidesFlag_func)
		{
			int i10 = 0, i11 = 0, i12 = 0, index13 = 0, index14 = 0, index15 = 0, index16 = 0;
			/*fix the position of y7 and y8*/

			i4 = 3;
			i12 = 7;
			/* try each permutation */
			for (i1 = 0; i1 < rowsRic_func ; i1++)
			if (i1 != i4 && i1 != i12)		
			for (i2 = 0; i2 < rowsRic_func; i2++)
			if (i2 != i1 && i2 != i4 && i2 != i12)						
			for (i3 = 0; i3 < rowsRic_func; i3++)
			if (i3 != i1 && i3 != i2 && i3 != i4 && i3 != i12) 
			for (i9 = 0; i9 < rowsRic_func; i9++)
			if (i9 != i1 && i9 != i2 && i9 != i3 && i9 != i4 && i9 != i12) 
			for (i10 = 0; i10 < rowsRic_func; i10++)
			if (i10 != i1 && i10 != i2 && i10 != i3 && i10 != i9 && i10 != i4 && i10 != i12) 
			for (i11 = 0; i11 < rowsRic_func; i11++)
			if (i11 != i1 && i11 != i2 && i11 != i3 && i11 != i9 && i11 != i10 && i11 != i4 && i11 != i12) 
			{
				if (permutationCounter > 5040)
					errorMessageAndExit("RiccatiMatrices_operations.c (find_permutation_with_minNorm)", 
										"permutation number higher than 720!!!!");
				
				index1 = i1;
				if(i1 >= rowsRic_func/2){index1 = i1 + rowsRic_func/2;}
				index2 = i2;
				if(i2 >= rowsRic_func/2){index2 = i2 + rowsRic_func/2;}
				index3 = i3;
				if(i3 >= rowsRic_func/2){index3 = i3 + rowsRic_func/2;}
				index4 = i4;
				
				index5 = index1 + rowsRic_func/2;
				index6 = index2 + rowsRic_func/2;
				index7 = index3 + rowsRic_func/2;
				index8 = index4 + rowsRic_func/2;
				
				index9 = i9;
				if(i9 >= rowsRic_func/2){index9 = i9 + rowsRic_func/2;}
				index10 = i10;
				if(i10 >= rowsRic_func/2){index10 = i10 + rowsRic_func/2;}
				index11 = i11;
				if(i11 >= rowsRic_func/2){index11 = i11 + rowsRic_func/2;}
				index12 = i12 + rowsRic_func/2;
				
				index13 = index9 + rowsRic_func/2;
				index14 = index10 + rowsRic_func/2;
				index15 = index11 + rowsRic_func/2;
				index16 = index12 + rowsRic_func/2;
				
				gsl_vector_set(vecDummySizeT_func,0,index1);
				gsl_vector_set(vecDummySizeT_func,1,index2);
				gsl_vector_set(vecDummySizeT_func,2,index3);
				gsl_vector_set(vecDummySizeT_func,3,index4);
				gsl_vector_set(vecDummySizeT_func,4,index5);
				gsl_vector_set(vecDummySizeT_func,5,index6);
				gsl_vector_set(vecDummySizeT_func,6,index7);
				gsl_vector_set(vecDummySizeT_func,7,index8);
				gsl_vector_set(vecDummySizeT_func,8,index9);
				gsl_vector_set(vecDummySizeT_func,9,index10);
				gsl_vector_set(vecDummySizeT_func,10,index11);
				gsl_vector_set(vecDummySizeT_func,11,index12);
				gsl_vector_set(vecDummySizeT_func,12,index13);
				gsl_vector_set(vecDummySizeT_func,13,index14);
				gsl_vector_set(vecDummySizeT_func,14,index15);
				gsl_vector_set(vecDummySizeT_func,15,index16);

				gsl_matrix_set_zero(matT_func);
				for (m = 0; m<rowsT_func; m++)
					gsl_matrix_set(matT_func,m, static_cast<int>(gsl_vector_get(vecDummySizeT_func,m)), 1.0);
				
				/********************************************************
				 ***** 
				 *****    performing check Eq (26) from Takata Loffler 2004
				 *****		     2004PASJ...56..645T
				 ***** 
				 *****	           det(T10*R + T11) != 0
				 *********************************************************/
				det_T10R_T11_s = calc_detRiccati_condition(matT_func, matT00_func, matT01_func, matT10_func,
														   matT11_func, matR_func, matDummySizeR_func, 
														   matDummySizeR_2_func, vecDummySizeR2_func, 
														   permDummySizeR_func, rowsRic_func);
				/*if the condition on the determinant is satisfied, then permute riccati*/
				if(det_T10R_T11_s != 0.0)
				{			
					permute_Riccati_matrix_R(matT_func, matR_func, matT00_func, matT01_func, matT10_func, matT11_func, 
											 matRperm_func, matDummySizeR_func, matDummySizeR_2_func, matDummySizeR_3_func, 
											 vecDummySizeR2_func, permDummySizeR_func, rowsRic_func);
					
					
					/* and now I have a permutation of R in matRperm_func. Calculate the norm of it, and 
					 store it together with the indices of the permutation. THis will be used later to find
					 the permutation with the minimum norm.*/						
					
					normR = 0.0;
					for (m=0; m<rowsRic_func; m++)
						for(n=0; n<rowsRic_func; n++){normR = normR + pow(gsl_matrix_get(matRperm_func,m,n), 2.0);}
					normR = sqrt(normR);
					
					gsl_matrix_set(mat_normR_func, k, 0,normR);
					
					/*store all permutations and corresponding ||R||*/
					for (m = 1; m<1+rowsT_func; m++)
						gsl_matrix_set(mat_normR_func,k, m, gsl_vector_get(vecDummySizeT_func,m-1));
					
					gsl_matrix_set(mat_normR_func, k, 1+rowsT_func,permutationCounter);
					
					k++;
				}//if(det_T10R_T11_s != 0.0)										
				permutationCounter++;
			}//if (i6 != i1 && i6 != i2 && i6 != i3 && i6 != i4 && i6 != i5)			
		}//if(tidesFlag_func)
		else
		{
			/* try each permutation */
			for (i1 = 0; i1 < rowsRic_func ; i1++)
			for (i2 = 0; i2 < rowsRic_func; i2++)
			if (i2 != i1)		
			for (i3 = 0; i3 < rowsRic_func; i3++)
			if (i3 != i1 && i3 != i2) 
			for (i7 = 0; i7 < rowsRic_func; i7++)
			if (i7 != i1 && i7 != i2 && i7 != i3)									
			for (i8 = 0; i8 < rowsRic_func; i8++)
			if (i8 != i1 && i8 != i2 && i8 != i3 && i8 != i7)									
			for (i9 = 0; i9 < rowsRic_func; i9++)
			if (i9 != i1 && i9 != i2 && i9 != i3 && i9 != i7 && i9 != i8)									
			{
				
				if (permutationCounter > 720)
					errorMessageAndExit("RiccatiMatrices_operations.c (find_permutation_with_minNorm)", 
										"permutation number higher than 720!!!!");
				
				index1 = i1;
				if(i1 >= rowsRic_func/2){index1 = i1 + rowsRic_func/2;}
				index2 = i2;
				if(i2 >= rowsRic_func/2){index2 = i2 + rowsRic_func/2;}
				index3 = i3;
				if(i3 >= rowsRic_func/2){index3 = i3 + rowsRic_func/2;}
				
				index4 = index1 + rowsRic_func/2;
				index5 = index2 + rowsRic_func/2;
				index6 = index3 + rowsRic_func/2;
				
				index7 = i7;
				if(i7 >= rowsRic_func/2){index7 = i7 + rowsRic_func/2;}
				index8 = i8;
				if(i8 >= rowsRic_func/2){index8 = i8 + rowsRic_func/2;}
				index9 = i9;
				if(i9 >= rowsRic_func/2){index9 = i9 + rowsRic_func/2;}
				
				index10 = index7 + rowsRic_func/2;
				index11 = index8 + rowsRic_func/2;
				index12 = index9 + rowsRic_func/2;
				
				gsl_vector_set(vecDummySizeT_func,0,index1);
				gsl_vector_set(vecDummySizeT_func,1,index2);
				gsl_vector_set(vecDummySizeT_func,2,index3);
				gsl_vector_set(vecDummySizeT_func,3,index4);
				gsl_vector_set(vecDummySizeT_func,4,index5);
				gsl_vector_set(vecDummySizeT_func,5,index6);
				gsl_vector_set(vecDummySizeT_func,6,index7);
				gsl_vector_set(vecDummySizeT_func,7,index8);
				gsl_vector_set(vecDummySizeT_func,8,index9);
				gsl_vector_set(vecDummySizeT_func,9,index10);
				gsl_vector_set(vecDummySizeT_func,10,index11);
				gsl_vector_set(vecDummySizeT_func,11,index12);
				
				gsl_matrix_set_zero(matT_func);
				for (m = 0; m<rowsT_func; m++)
					gsl_matrix_set(matT_func,m, static_cast<int>(gsl_vector_get(vecDummySizeT_func,m)), 1.0);
				
				/********************************************************
				 ***** 
				 *****    performing check Eq (26) from Takata Loffler 2004
				 *****		     2004PASJ...56..645T
				 ***** 
				 *****	           det(T10*R + T11) != 0
				 *********************************************************/
				det_T10R_T11_s = calc_detRiccati_condition(matT_func, matT00_func, matT01_func, matT10_func,
														   matT11_func, matR_func, matDummySizeR_func, 
														   matDummySizeR_2_func, vecDummySizeR2_func, 
														   permDummySizeR_func, rowsRic_func);
				
				/*if the condition on the determinant is satisfied, then permute riccati*/
				if(det_T10R_T11_s != 0.0)
				{			
					
					permute_Riccati_matrix_R(matT_func, matR_func, matT00_func, matT01_func, matT10_func, matT11_func, 
											 matRperm_func, matDummySizeR_func, matDummySizeR_2_func, matDummySizeR_3_func, 
											 vecDummySizeR2_func, permDummySizeR_func, rowsRic_func);
					
					
					/* and now I have a permutation of R in matRperm_func. Calculate the norm of it, and 
					 store it together with the indices of the permutation. THis will be used later to find
					 the permutation with the minimum norm.*/						
					
					normR = 0.0;
					for (m=0; m<rowsRic_func; m++)
						for(n=0; n<rowsRic_func; n++){normR = normR + pow(gsl_matrix_get(matRperm_func,m,n), 2.0);}
					normR = sqrt(normR);
					
					gsl_matrix_set(mat_normR_func, k, 0,normR);
					
					/*store all permutations and corresponding ||R||*/
					for (m = 1; m<1+rowsT_func; m++)
						gsl_matrix_set(mat_normR_func,k, m, gsl_vector_get(vecDummySizeT_func,m-1));
					
					gsl_matrix_set(mat_normR_func, k, 1+rowsT_func,permutationCounter);
					
					k++;
				}//if(det_T10R_T11_s != 0.0)										
				permutationCounter++;
			}//if (i6 != i1 && i6 != i2 && i6 != i3 && i6 != i4 && i6 != i5)
			
		}//else non tides	
		
	}//	if non adiabatic
	else
		errorMessageAndExit("RiccatiMatrices_operations.c (find_permutation_with_minNorm)", 
							"T with minimum ||R|| not calculated. Adiabatic or what?!");
		
	/*number of elements that satisfied the condition on the determinant*/	
	normR_counts = k;										
	
	/*find element with smallest norm, and apply that trasformation*/
	normR_min = gsl_matrix_get(mat_normR_func, 0, 0);

	/*assuming that the first element in mat_normR_func is the one with minimum normR, 
	 and storing the indices in VecIndices*/
	for (m = 0; m<rowsT_func; m++)
		gsl_vector_set(vecDummySizeT_func,m,gsl_matrix_get(mat_normR_func, 0, m+1));
	
	/*Looping on all permutations applied*/
	for (k = 1; k<normR_counts; k++)
	{		
		/*if a lower norm is found,the corresponding indices are stored in VecIndices*/
		if (gsl_matrix_get(mat_normR_func, k, 0) <= normR_min)
		{
			normR_min = gsl_matrix_get(mat_normR_func, k, 0);
			PermutationMinNormNumber_func = gsl_matrix_get(mat_normR_func, k, 1+rowsT_func);
			
			for (m = 0; m<rowsT_func; m++)
				gsl_vector_set(vecDummySizeT_func,m,gsl_matrix_get(mat_normR_func, k, m+1));
		}	
	}
	
	fill_permutation_T(matT_func, vecDummySizeT_func, rowsT_func);	
	
	gsl_vector_memcpy (vecPermutIndices_func, vecDummySizeT_func);
	
	free(mat_normR_func);
	free(matRperm_func);
}//find_permutation_with_minNorm
/****************************************************************************************
 ****************************************************************************************/
void fill_permutation_T(gsl_matrix *matT_func, gsl_vector *vecPermutIndices_func, int rowsT_func)
{
	int m = 0;
	gsl_matrix_set_zero(matT_func);		
	
	for (m = 0; m<rowsT_func; m++)
		gsl_matrix_set(matT_func,m, static_cast<int>(gsl_vector_get(vecPermutIndices_func,m)), 1.0);
	
	
}//fill_permutation_T


void permute_Riccati_UV(gsl_vector *vecU_func, gsl_vector *vecV_func, gsl_matrix *Rperm_func, 
						gsl_matrix *T_func, gsl_matrix *T00_func, gsl_matrix *T01_func, 
						gsl_matrix *T10_func, gsl_matrix *T11_func, gsl_vector *vecDummySizeR_1_func, 
						gsl_vector *vecDummySizeR_2_func, int rowsRic_func)
{
	
	/********************************************************
	 ***** 
	 *****    Applying Eq (24) from Takata Loffler 2004
	 *****		     2004PASJ...56..645T
	 *****		V' = T10U + T11V
	 *****		and calculating U' from V' and R'
	 *****     
	 *********************************************************/
	extract_submatrices(T_func, T00_func, T01_func, T10_func, T11_func, rowsRic_func);

	/* T10U --> vecDummySizeR_1_func */
	mul_matrix_by_vec(T10_func, vecU_func, vecDummySizeR_1_func, vecDummySizeR_1_func, rowsRic_func);

	/* T11V --> vecDummySizeR_2_func */
	mul_matrix_by_vec(T11_func, vecV_func, vecDummySizeR_2_func, vecDummySizeR_2_func, rowsRic_func);

	/* creating V' */
	vector_sum (vecDummySizeR_1_func, vecDummySizeR_2_func, vecV_func, rowsRic_func);

	/* creating U' = R'V' */
	mul_matrix_by_vec(Rperm_func, vecV_func, vecU_func, vecU_func, rowsRic_func);
}//permute_Riccati_UV


void back_to_original_Riccati_UV(gsl_matrix *T_func, gsl_matrix *T00_func, gsl_matrix *T01_func, 
								 gsl_matrix *T10_func, gsl_matrix *T11_func, 
								 gsl_vector *U_func, gsl_vector *V_func, gsl_vector *Uoriginal_func, 
								 gsl_vector *Voriginal_func, gsl_matrix *matDummySizeT_func, 
								 gsl_matrix *matDummySizeT_2_func, gsl_vector *vecDummySizeR_1_func, 
								 gsl_vector *vecDummySizeR_2_func, gsl_vector *vecDummySizeT_func,
								 gsl_permutation *permDummySizeT_func, int rowsRic_func)
{
	
	/********************************************************
	 ***** 
	 *****    Inverting Eq (24) from Takata Loffler 2004
	 *****		     2004PASJ...56..645T
	 ***** 
	 *****	     T^-1(U',V') = (U,V)
	 *****	     
	 *****	     U = T00U'+T01V'
	 *****	     V = T10U'+T11V'
	 *****	     
	 *****	     With Tij submatrices of the inverted
	 *****	     matrix T.
	 *********************************************************/
	int i = 0, rowsT_func = 2*rowsRic_func;

	/* Calculate T^-1 */
	InverseMatrix(T_func, matDummySizeT_func, rowsT_func, permDummySizeT_func, matDummySizeT_2_func);
	
	/* extract submatrices from inverted T --> T00, T01, T10, T11 */
	extract_submatrices(matDummySizeT_func, T00_func, T01_func, T10_func, T11_func, rowsRic_func);
	
	gsl_vector_set_zero(vecDummySizeT_func);
	
	/* T00 U' --> vecDummysizeR_1_func */
	mul_matrix_by_vec(T00_func, U_func, vecDummySizeR_1_func,vecDummySizeR_1_func, rowsRic_func);

	/* T01 V' --> vecDummysizeR_2_func */
	mul_matrix_by_vec(T01_func, V_func, vecDummySizeR_2_func, vecDummySizeR_2_func, rowsRic_func);

	/* T00 U' + T01 V' in a vector */
	for (i = 0; i<rowsRic_func; i++)
		gsl_vector_set(vecDummySizeT_func,i,
					   gsl_vector_get(vecDummySizeR_1_func,i)+gsl_vector_get(vecDummySizeR_2_func,i));	

	/* T10 U' --> vecDummysizeR_1_func */
	mul_matrix_by_vec(T10_func, U_func, vecDummySizeR_1_func, vecDummySizeR_1_func, rowsRic_func);

	/* T11 V' --> vecDummysizeR_1_func */
	mul_matrix_by_vec(T11_func, V_func, vecDummySizeR_2_func, vecDummySizeR_2_func, rowsRic_func);

	/* T10 U'+T11 V' in a vector */
	for (i = 0; i<rowsRic_func; i++)
		gsl_vector_set(vecDummySizeT_func,i+rowsRic_func,
					   gsl_vector_get(vecDummySizeR_1_func,i)+gsl_vector_get(vecDummySizeR_2_func,i));

	/* create U and V original */
	for (i = 0; i<rowsRic_func; i++)
	{
		gsl_vector_set(Uoriginal_func,i, gsl_vector_get(vecDummySizeT_func,i));
		gsl_vector_set(Voriginal_func,i, gsl_vector_get(vecDummySizeT_func,i+rowsRic_func));
	}
}//back_to_original_Riccati_UV

void permute_Riccati_V(gsl_vector *vecU_func, gsl_vector *vecV_func, gsl_matrix *T_func, 
						gsl_matrix *T00_func, gsl_matrix *T01_func, gsl_matrix *T10_func, gsl_matrix *T11_func, 
						gsl_vector *vecDummySizeR_1_func, gsl_vector *vecDummySizeR_2_func, int rowsRic_func)
{
	
	/********************************************************
	 ***** 
	 *****    Applying Eq (24) from Takata Loffler 2004
	 *****		     2004PASJ...56..645T
	 *****		v' = T10U + T11V
	 *********************************************************/
	extract_submatrices(T_func, T00_func, T01_func, T10_func, T11_func, rowsRic_func);
	
	/* T10 U --> vecDummysizeR_1_func */
	mul_matrix_by_vec(T10_func, vecU_func, vecDummySizeR_1_func, vecDummySizeR_1_func, rowsRic_func);

	/* T11 V --> vecDummysizeR_2_func */
	mul_matrix_by_vec(T11_func, vecV_func, vecDummySizeR_2_func, vecDummySizeR_2_func, rowsRic_func);
	
	/* Create V' */
	vector_sum (vecDummySizeR_1_func, vecDummySizeR_2_func, vecV_func, rowsRic_func);	
}//permute_Riccati_V

void store_IC_ForIntegrationEigenfunctions(gsl_matrix *matPermInfo_func, double csi_i_func, double csi_f_func, 
										   gsl_vector *permutInit_BC_tot_func, 
										   gsl_vector *rijPrime_func, 
										   gsl_vector *permutInit_BC_step_func,  
										   int &counterPermut_func, int sizeR_func)
{
	/*Store initial conditions in matrix that will be used during integration of
	 eigenfunctions. See beginning of the main file for an explanation of what 
	 matPermInfo_func contains*/
	 
	int rowsT_func= static_cast<int>(sqrt(sizeR_func)*2), 
	i = 0, m = 2, locationOf1 = 0;
	double r_ij = 0.0;

	counterPermut_func = 0;
	gsl_matrix_set_zero(matPermInfo_func);
	gsl_matrix_set(matPermInfo_func,counterPermut_func, 0, csi_i_func);
	gsl_matrix_set(matPermInfo_func,counterPermut_func, 1, csi_f_func);
	
	/*Store initial permutations total and partial*/
	for(i = 0; i<rowsT_func; i++)
	{
		locationOf1 = static_cast<int>(gsl_vector_get(permutInit_BC_tot_func,i));
		gsl_matrix_set(matPermInfo_func,counterPermut_func, i+m, locationOf1);
		locationOf1 = static_cast<int>(gsl_vector_get(permutInit_BC_step_func,i));
		gsl_matrix_set(matPermInfo_func,counterPermut_func, i+m+rowsT_func+sizeR_func, locationOf1);
	}
	
	/*Store riccati components.*/
	for(i = 0; i<sizeR_func; i++)
	{
		r_ij = gsl_vector_get(rijPrime_func,i);
		gsl_matrix_set(matPermInfo_func,counterPermut_func, i+m+rowsT_func, r_ij);	
	}
	
}//store_IC_ForIntegrationEigenfunctions

void update_csi_and_r_ij_forEndInterval(gsl_matrix *matPermInfo_func, double csi_f_func, 
										gsl_vector *rijPrime_func, int counterPermut_func, int sizeR_func)
{
	int rowsT_func= static_cast<int>(sqrt(sizeR_func)*2), 
	i = 0, i_start = 0, i_end = 0;
	double r_ij = 0.0;

	/*counting the first two elements (integration extremes for sub-interval), and 
	 the indices marking the position of "1s" in Ttot*/
	i_start = 2 + rowsT_func;
	i_end = i_start + sizeR_func;
	
	gsl_matrix_set(matPermInfo_func,counterPermut_func, 1, csi_f_func);
	
	for (i = i_start; i<i_end; i++)
	{
		r_ij = gsl_vector_get(rijPrime_func, i - i_start);
		gsl_matrix_set(matPermInfo_func,counterPermut_func, i, r_ij);
	}
	
}//update_csi_and_r_ij_forEndInterval

void store_Tstep_ForIntegrationEigenfunctions(gsl_matrix *matPermInfo_func, 
											  gsl_vector *vecPermutIndices_func,
											  int sizeR_func, 
											  int &counterPermut_func, int PermutationMinNormNumber_func)
{
	
	int rowsT_func= static_cast<int>(sqrt(sizeR_func)*2), 
	i = 0, 
	i_start = 2+rowsT_func + sizeR_func, 
	i_end = 2+rowsT_func + sizeR_func + rowsT_func;
	counterPermut_func++;
	
	for (i = i_start; i< i_end; i++)
		gsl_matrix_set(matPermInfo_func, counterPermut_func, i, gsl_vector_get(vecPermutIndices_func,i - i_start));
	
	gsl_matrix_set(matPermInfo_func, counterPermut_func, i_end, PermutationMinNormNumber_func);
	
}//store_Tstep_ForIntegrationEigenfunctions


void store_Ttot_ForIntegrationEigenfunctions(gsl_matrix *matT_func, 
											 int rowsT_func, int colsT_func, 
											 gsl_matrix *matPermInfo_func, 
											 int counterPermut_func)
{
	int t_ij = 0, i = 0, j = 0;
	
	for (i=0; i<rowsT_func; i++)
		for (j=0; j<colsT_func; j++)
		{
			t_ij = static_cast<int>(gsl_matrix_get(matT_func,i,j));
			if (t_ij)
				gsl_matrix_set(matPermInfo_func, counterPermut_func, 2+i, j);
		}
}//store_Ttot_ForIntegrationEigenfunctions



/*************************************************************************************
 ************ 
 ************	These routines creates the permuted C and T for each possible
 ************	permutation matrix T (24 in total)
 ************ 
 **************************************************************************************/

void create_permutedCD(double csi_func, double Vg_func,double Astar_func, double U_func,double l_func,						 
					   double omega_r_func, double omega_i_func, double Vt_func, double V_func, double del_func, double delad_func, 							
					   double ks_func, double c1_func, double c2_func, double c3_func, double c4_func, double epsAd_func, double dlnLR_dlnr_func, 
					   double epsS_func, gsl_matrix *C_func, gsl_matrix *D_func, int rowsRic_func, int permutationNumber_func)
{
	gsl_matrix_set_zero(C_func);
	gsl_matrix_set_zero(D_func);
	
	int i = 0, j = 0;
	double llP1 = l_func*(1.0 + l_func), w_r2 = omega_r_func*omega_r_func, w_i2 = omega_i_func*omega_i_func;

	if(permutationNumber_func == 569)
	{
		gsl_matrix_set(C_func,0, 0, 0.0);
		gsl_matrix_set(C_func,0, 1, (delad_func*llP1)/del_func - c3_func*epsAd_func*V_func);
		gsl_matrix_set(C_func,0, 2, (-llP1 + c3_func*del_func*epsS_func*V_func + c4_func*del_func*V_func*omega_i_func)/(del_func*V_func));		
		gsl_matrix_set(C_func,0, 3, 0.0);
		gsl_matrix_set(C_func,0, 4, 0.0);
		gsl_matrix_set(C_func,0, 5, c4_func*omega_r_func);
		gsl_matrix_set(C_func,1, 0, 0.0);
		gsl_matrix_set(C_func,1, 1, -Astar_func);
		gsl_matrix_set(C_func,1, 2, Vt_func);
		gsl_matrix_set(C_func,1, 3, 0.0);
		gsl_matrix_set(C_func,1, 4, 0.0);
		gsl_matrix_set(C_func,1, 5, 0.0);
		gsl_matrix_set(C_func,2, 0, 0.0);
		gsl_matrix_set(C_func,2, 1, Vg_func);
		gsl_matrix_set(C_func,2, 2, Vt_func);		
		gsl_matrix_set(C_func,2, 3, 0.0);
		gsl_matrix_set(C_func,2, 4, 0.0);
		gsl_matrix_set(C_func,2, 5, 0.0);		
		
		gsl_matrix_set(D_func,0, 0, 2.0 - dlnLR_dlnr_func - l_func);
		gsl_matrix_set(D_func,0, 1, (-(delad_func*llP1) + c3_func*del_func*epsAd_func*V_func)/del_func - 
					   (c3_func*llP1*(w_i2 - w_r2))/(c1_func*(w_i2 + w_r2)*(w_i2 + w_r2)));
		gsl_matrix_set(D_func,0, 2, -(((-delad_func + del_func)*llP1)/del_func) - c3_func*epsAd_func*V_func);
		gsl_matrix_set(D_func,0, 3, 0);
		gsl_matrix_set(D_func,0, 4, (2*c3_func*llP1*omega_i_func*omega_r_func)/(c1_func*(w_i2 + w_r2)*(w_i2 + w_r2)));
		gsl_matrix_set(D_func,0, 5, 0);
		gsl_matrix_set(D_func,1, 0, 0);
		gsl_matrix_set(D_func,1, 1, 3.0 + Astar_func - l_func - U_func);
		gsl_matrix_set(D_func,1, 2, -Astar_func + c1_func*(-w_i2 + w_r2));
		gsl_matrix_set(D_func,1, 3, 0);
		gsl_matrix_set(D_func,1, 4, 0);
		gsl_matrix_set(D_func,1, 5, -2.0*c1_func*omega_i_func*omega_r_func);
		gsl_matrix_set(D_func,2, 0, 0);
		gsl_matrix_set(D_func,2, 1, -Vg_func + (llP1*(-w_i2 + w_r2))/(c1_func*(w_i2 + w_r2)*(w_i2 + w_r2)));
		gsl_matrix_set(D_func,2, 2, -1.0 - l_func + Vg_func);
		gsl_matrix_set(D_func,2, 3, 0);
		gsl_matrix_set(D_func,2, 4, (2.0*llP1*omega_i_func*omega_r_func)/(c1_func*(w_i2 + w_r2)*(w_i2 + w_r2)));
		gsl_matrix_set(D_func,2, 5, 0);				
	}//if(PermutationMinNormNum[[2ber_func == 569)
	else if(permutationNumber_func == 1000)
	{
		gsl_matrix_set(C_func,0, 0, -Astar_func + c1_func*(-w_i2 + w_r2));
		gsl_matrix_set(C_func,0, 1, 0);
		gsl_matrix_set(C_func,0, 2, Vt_func);
		gsl_matrix_set(C_func,0, 3, -2.0*c1_func*omega_i_func*omega_r_func);
		gsl_matrix_set(C_func,0, 4, 0);
		gsl_matrix_set(C_func,0, 5, 0);		
		gsl_matrix_set(C_func,1, 0, 0);		
		gsl_matrix_set(C_func,1, 1, 1);		
		gsl_matrix_set(C_func,1, 2, 0);		
		gsl_matrix_set(C_func,1, 3, 0);		
		gsl_matrix_set(C_func,1, 4, 0);		
		gsl_matrix_set(C_func,1, 5, 0);		
		gsl_matrix_set(C_func,2, 0, -(((-delad_func + del_func)*llP1)/del_func) - 
					   c3_func*epsAd_func*V_func);
		gsl_matrix_set(C_func,2, 1, 0);
		gsl_matrix_set(C_func,2, 2, (-llP1 + c3_func*del_func*epsS_func*V_func + c4_func*del_func*V_func*omega_i_func)/(del_func*V_func));
		gsl_matrix_set(C_func,2, 3, 0);
		gsl_matrix_set(C_func,2, 4, 0);
		gsl_matrix_set(C_func,2, 5, c4_func*omega_r_func);
		
		gsl_matrix_set(D_func,0, 0, 3.0 + Astar_func - l_func - U_func);
		gsl_matrix_set(D_func,0, 1, -Astar_func);
		gsl_matrix_set(D_func,0, 2, 0);
		gsl_matrix_set(D_func,0, 3, 0);
		gsl_matrix_set(D_func,0, 4, 0);
		gsl_matrix_set(D_func,0, 5, 0);
		gsl_matrix_set(D_func,1, 0, 0);
		gsl_matrix_set(D_func,1, 1, 3.0 - l_func - U_func);
		gsl_matrix_set(D_func,1, 2, 0);
		gsl_matrix_set(D_func,1, 3, 0);
		gsl_matrix_set(D_func,1, 4, 0);
		gsl_matrix_set(D_func,1, 5, 0);
		gsl_matrix_set(D_func,2, 0, (-(delad_func*llP1) + c3_func*del_func*epsAd_func*V_func)/del_func - 
					   (c3_func*llP1*(w_i2 - w_r2))/(c1_func*(w_i2 + w_r2)*(w_i2 + w_r2)));
		gsl_matrix_set(D_func,2, 1, (delad_func*llP1)/del_func - c3_func*epsAd_func*V_func);
		gsl_matrix_set(D_func,2, 2, 2.0 - dlnLR_dlnr_func - l_func);
		gsl_matrix_set(D_func,2, 3, (2.0*c3_func*llP1*omega_i_func*omega_r_func)/(c1_func*(w_i2 + w_r2)*(w_i2 + w_r2)));
		gsl_matrix_set(D_func,2, 4, 0);
		gsl_matrix_set(D_func,2, 5, 0);
	}//else if(PermutationMinNormNumber_func == 1000)
	else
		errorMessageAndExit("RiccatiMatrices_operations.c (create_permutedCD)", "unkown permutationNumber_func!");
		
		gsl_matrix_set(C_func,3, 0, -gsl_matrix_get(C_func,0, 3));
		gsl_matrix_set(C_func,3, 1, -gsl_matrix_get(C_func,0, 4));
		gsl_matrix_set(C_func,3, 2, -gsl_matrix_get(C_func,0, 5));
		gsl_matrix_set(C_func,3, 3, gsl_matrix_get(C_func,0, 0));
		gsl_matrix_set(C_func,3, 4, gsl_matrix_get(C_func,0, 1));
		gsl_matrix_set(C_func,3, 5, gsl_matrix_get(C_func,0, 2));		
		gsl_matrix_set(C_func,4, 0, -gsl_matrix_get(C_func,1, 3));
		gsl_matrix_set(C_func,4, 1, -gsl_matrix_get(C_func,1, 4));
		gsl_matrix_set(C_func,4, 2, -gsl_matrix_get(C_func,1, 5));
		gsl_matrix_set(C_func,4, 3, gsl_matrix_get(C_func,1, 0));
		gsl_matrix_set(C_func,4, 4, gsl_matrix_get(C_func,1, 1));
		gsl_matrix_set(C_func,4, 5, gsl_matrix_get(C_func,1, 2));		
		gsl_matrix_set(C_func,5, 0, -gsl_matrix_get(C_func,2, 3));
		gsl_matrix_set(C_func,5, 1, -gsl_matrix_get(C_func,2, 4));
		gsl_matrix_set(C_func,5, 2, -gsl_matrix_get(C_func,2, 5));
		gsl_matrix_set(C_func,5, 3, gsl_matrix_get(C_func,2, 0));
		gsl_matrix_set(C_func,5, 4, gsl_matrix_get(C_func,2, 1));
		gsl_matrix_set(C_func,5, 5, gsl_matrix_get(C_func,2, 2));
		
		gsl_matrix_set(D_func,3, 0, -gsl_matrix_get(D_func,0, 3));
		gsl_matrix_set(D_func,3, 1, -gsl_matrix_get(D_func,0, 4));
		gsl_matrix_set(D_func,3, 2, -gsl_matrix_get(D_func,0, 5));
		gsl_matrix_set(D_func,3, 3, gsl_matrix_get(D_func,0, 0));
		gsl_matrix_set(D_func,3, 4, gsl_matrix_get(D_func,0, 1));
		gsl_matrix_set(D_func,3, 5, gsl_matrix_get(D_func,0, 2));
		gsl_matrix_set(D_func,4, 0, -gsl_matrix_get(D_func,1, 3));
		gsl_matrix_set(D_func,4, 1, -gsl_matrix_get(D_func,1, 4));
		gsl_matrix_set(D_func,4, 2, -gsl_matrix_get(D_func,1, 5));
		gsl_matrix_set(D_func,4, 3, gsl_matrix_get(D_func,1, 0));
		gsl_matrix_set(D_func,4, 4, gsl_matrix_get(D_func,1, 1));
		gsl_matrix_set(D_func,4, 5, gsl_matrix_get(D_func,1, 2));		
		gsl_matrix_set(D_func,5, 0, -gsl_matrix_get(D_func,2, 3));
		gsl_matrix_set(D_func,5, 1, -gsl_matrix_get(D_func,2, 4));
		gsl_matrix_set(D_func,5, 2, -gsl_matrix_get(D_func,2, 5));
		gsl_matrix_set(D_func,5, 3, gsl_matrix_get(D_func,2, 0));
		gsl_matrix_set(D_func,5, 4, gsl_matrix_get(D_func,2, 1));
		gsl_matrix_set(D_func,5, 5, gsl_matrix_get(D_func,2, 2));
		
		for (i = 0; i<rowsRic_func; i++)
			for(j = 0; j<rowsRic_func; j++)
			{
				gsl_matrix_set(C_func, i, j, gsl_matrix_get(C_func, i, j)/csi_func);
				gsl_matrix_set(D_func, i, j, gsl_matrix_get(D_func, i, j)/csi_func);
				
			}	
}//create_permutedCD

/*************************************************************************************
 ************	Complete R using the properties of simmetry of R for 
 ************   the non adiabatic case.
 ************   
 ************  
 ************   R = 
 ************   |R00   R01|    |R00   R01|
 ************   |R10   R11| =  |-R01  R00| 
 ************  
 ************   for instance, if R = (6 x 6) 
 ************   |0    1    2    3    4    5|    
 ************   |6    7    8    9   10   11|
 ************   |12  13   14   15   16   17|
 ************   |18  19   20   21   22   23| 
 ************   |24  25   26   27   28   29|
 ************   |30  31   32   33   34   35|
 ************  
 ************   
 **************************************************************************************/

void complete_y_Riccati_usingSimmetryR(double *y_Riccati_func, int rowsR_func)
{
	
	int k = 0, i = 0, i_in = 0, offset = 0;
	
	/*using simmetry between R01 and R10*/
	
	/*The offset is the starting point for filling. In the R = (6 x 6) above offset it 18*/
	offset = rowsR_func*rowsR_func/2;	
	for(k = 0; k < rowsR_func/2; k++)   
	{
		/*determining the row position for R = (6 x 6) above i_in is 18, 24, 30*/
		i_in =  offset + k * rowsR_func;		
		
		/*Filling in each row half the way through. for R = (6 x 6) this fills in elements (18, 19, 20), (24, 25, 26), (30, 31, 32)
		 using the elements (3, 4, 5), (9, 10, 11), (15, 16, 17)*/
		
		for(i = i_in; i < i_in + rowsR_func/2; i++)  
			y_Riccati_func[i] = -y_Riccati_func[i - offset+rowsR_func/2];
	}
	
	/*using simmetry between R00 and R11*/
	
	/* In the R = (6 x 6) case above, offset is moved by 3 more to reach the down-right sub matrix --> 21*/
	
	offset = rowsR_func*rowsR_func/2 + rowsR_func/2;
	for(k = 0; k < rowsR_func/2; k++)
	{
		i_in =  offset + k * rowsR_func;
		for(i = i_in; i < i_in + rowsR_func/2; i++)
			y_Riccati_func[i] = y_Riccati_func[i - offset];
	}
	
	
}//void complete_y_Riccati_usingSimmetryR
/*************************************************************************************
*************************************************************************************/
void complete_matR_usingSymmetry(gsl_matrix *matR_func, int rowsR_func)
{
	int halfRowR = static_cast<int>(rowsR_func/2), i = 0, j = 0;
	
	for(i = halfRowR; i<rowsR_func; i++)
		for(j = halfRowR; j<rowsR_func; j++)
		{
			gsl_matrix_set(matR_func,i,j, gsl_matrix_get(matR_func, i-halfRowR, j-halfRowR));		
			
			gsl_matrix_set(matR_func,i,j-halfRowR, -gsl_matrix_get(matR_func, i-halfRowR, j));	
		}
}//void complete_matR_usingSymmetry



