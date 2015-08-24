/**************************************************************************
 **************************************************************************
 ************           
 ************    Integrator for eigenfunctions and eigenfrequencies
 ************    following Takata & Loffler 2004
 ************              2004PASJ...56..645T
 ************           
 ************    See Eq.(13) and Eq. (22)
 ************     
 ************     dR/dx = B+AR-RD-RCR
 ************     dv/dx = (CR+D)v
 ************     
 ************     
 ************     Variables:
 ************     
 **************************************************************************
 ************     FOr adiabatic: func_Riccati_R and jac_Riccati_R are used.
 ************     y[0] = r00
 ************     y[1] = r01
 ************     y[2] = r10
 ************     y[3] = r11
 ************     y[4] = v0
 ************     y[5] = v1	 
 ************     
 ************     Procedure:
 ************     --Integration of R to find the EIGENFREQUENCIES--
 ************     
 ************     Rij are integrated, and Vi are sobstituted with the first row of the matrix R (better than setting
 ************     the derivatives to zero
 ************     
 ************     --Integration of V to find the EIGENFUNCTIONS--
 ************     I have two options:
 ************     1) Integrating both R and V at the same time, which could be unstable at the star's boundaries
 ************     2) Integrate V and interpolate R.In this sencond case, the derivatives of R are
 ************     substituted with the derivatives of V
 ************     (i.e. dR00 = dV0, dR01 = dV1
 ************     dR10 = dV0, dR11 = dV1
 ************     dR20 = dV0, dR21 = dV1)
 ************     and the Jacobian is re-calculated accordingly.
 ************     
 ************     
 **************************************************************************
 ************     FOr non-adiabatic: func_Riccati_Ronly, jac_Riccati_Ronly,
 ************     func_Riccati_Vonly and jac_Riccati_Vonly are used
 ************     
 ************     
 ************     Func_Riccati_Ronly:
 ************     y[0] = r00	y[6] = r10	y[12] = r20    y[18] = r30	y[24] = r40	y[30] = r50     
 ************     y[1] = r01	y[7] = r11	y[13] = r21    y[19] = r31	y[25] = r41	y[31] = r51     
 ************     y[2] = r02	y[8] = r12	y[14] = r22    y[20] = r32	y[26] = r42	y[32] = r52     
 ************     y[3] = r03	y[9] = r13	y[15] = r23    y[21] = r33	y[27] = r43	y[33] = r53     
 ************     y[4] = r04	y[10] = r14	y[16] = r24    y[22] = r34	y[28] = r44	y[34] = r54     
 ************     y[5] = r05	y[11] = r15	y[17] = r25    y[23] = r35	y[29] = r45	y[35] = r55     
 ************     
 ************     
 ************     Func_Riccati_Vonly:
 ************     y[0] = v0	
 ************     y[1] = v1	
 ************     y[2] = v2	
 ************     y[3] = v3	
 ************     y[4] = v4	
 ************     y[5] = v5	
 ************     
 ************     Procedure:
 ************     --Integration of R to find the EIGENFREQUENCIES--
 ************     
 ************     Rij are integrated
 ************     
 ************     --Integration of V to find the EIGENFUNCTIONS--
 ************     
 ************     Vi are integrated while Rij are interpolated
 ************     
 **************************************************************************
 *************************************************************************/

#include <iomanip>
#include <vector>
#include <iostream>
#include <math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "gsl_odeiv_evolve_Riccati.h"
#include "params_integrator_V_nonAd_rkf45.h"
#include "params_integrator_V_nonAd.h"
#include "params_integrator_adiabatic.h"
#include "params_integrator_R_nonAd.h"
#include "RiccatiMatrices_operations.h"
#include "SteffenInterp.h"
#include "matrix_operations.h"
#include "eigenfunction_operations.h"
#include "IOfiles.h"

using namespace std;


/********************************************************************************
 ********************************************************************************
 ******************    
 ******************				Adiabatic case:
 ******************    
 ********************************************************************************
 ********************************************************************************/					

int func_Riccati_R (double csi, const double y[], double f[], void *params)
{	
	
	int i = 0, j = 0, k = 0, m = 0;
  	params_integrator_adiabatic_struct *p = (params_integrator_adiabatic_struct *) params;
	double omega_r 	 = p->omega_r_struct,
	l   	         = p->l_struct;
	int		polyMesh = p->polyMesh_struct,
	rowsRic	         = p->rowsRic_struct,
	interpolate_rij  = p->interpol_rij_struct, 
	sizeR = rowsRic*rowsRic; 
	const int position_R_y_Riccati = 0, 
	position_V_y_Riccati = sizeR;

	double lowEndCsiCalc = 0.0, highEndCsiCalc = 0.0;
	int integratingV = 0, integratingV_fitSurface = 0, integratingV_fitCenter = 0; 
	
	/********************************************************************************
	 ***    Determining what I am integrating
	 ********************************************************************************/					
	if(p->inWard_outWard_struct <= 1){integratingV = 1;}
	if(p->inWard_outWard_struct == 0){integratingV_fitSurface = 1;}
	if(p->inWard_outWard_struct == 1){integratingV_fitCenter = 1;}
	
	/********************************************************************************
	 ***    forcing csi to be within integration interval
	 ***    to avoid numerical errors to the 16th decimal place
	 ********************************************************************************/				
	if (integratingV_fitSurface)
	{
		if (csi > p->csiFin_force_struct){csi = p->csiFin_force_struct;}
		if (csi < p->csiIn_force_struct){csi = p->csiIn_force_struct;}
	}
	if (integratingV_fitCenter)
	{
		if (csi < p->csiFin_force_struct){csi = p->csiFin_force_struct;}
		if (csi > p->csiIn_force_struct){csi = p->csiIn_force_struct;}
	}
	/********************************************************************************
	 ***    Interpolating with Steffen the stellar model's quantities
	 ********************************************************************************/				
	evalSteffenInterp(polyMesh,p->csiRel_vec_struct,p->fitVgCoeffs_struct,csi, p->VgFuncs_struct);	
	evalSteffenInterp(polyMesh,p->csiRel_vec_struct,p->fitAstarCoeffs_struct,csi, p->AstarFuncs_struct);
	evalSteffenInterp(polyMesh,p->csiRel_vec_struct,p->fitUcoeffs_struct,csi, p->UFuncs_struct);
	evalSteffenInterp(polyMesh,p->csiRel_vec_struct,p->fitC1coeffs_struct,csi, p->c1Funcs_struct); 

	/********************************************************************************
	 ***    Starting from original ABCD, going to M, permuting M and back to ABCD
	 ***    after this:T and M are unchanged, while A, B, C, D are permuted 
	 ********************************************************************************/				
	permute_Riccati_ABCD(csi, p->VgFuncs_struct[0], p->AstarFuncs_struct[0], 
						 p->UFuncs_struct[0], p->c1Funcs_struct[0], l, omega_r, 						 
						 p->omega_i_struct, p->VtFuncs_struct[0], p->VFuncs_struct[0], 
						 p->delFuncs_struct[0], p->delADFuncs_struct[0], p->ksFuncs_struct[0], 
						 p->c2Funcs_struct[0], p->c3Funcs_struct[0], p->c4Funcs_struct[0], 
						 p->epsADFuncs_struct[0], p->dlnLR_dlnrFuncs_struct[0], p->epsSFuncs_struct[0], 
						 p->matA_struct, p->matB_struct, p->matC_struct,p->matD_struct, 
						 p->matM_struct, p->matT_struct,
						 p-> matDummySizeT_struct, p-> matDummySizeT_2_struct, 
						 p->vecDummySizeT2_struct, p->permDummySizeT_struct, rowsRic, 
						 p->inWard_outWard_struct, p->nonAdiabatic_struct, p->tidesFlag_struct);	
	/********************************************************************************
	 ***  Fill in Riccati matrix to create Riccati equations ==>
	 ***  Matrix R contains the calculated values
	 ********************************************************************************/					
	m=0;
	for (k=0; k<rowsRic; k++)
		for (j=0; j<rowsRic; j++, m++)
			gsl_matrix_set(p->matR_struct,k,j, y[m+position_R_y_Riccati]);	

	/********************************************************************************
	 ***    Create Riccati equations:
	 ***    after this:A, B, C, D and R are unchanged, while i've dR/dx and
	 ***    R still contains the calculated values
	 ********************************************************************************/				
	create_Riccati_eqs(p -> matA_struct, p -> matB_struct, p -> matC_struct, p -> matD_struct, 
					   p -> matR_struct, p -> matdR_struct, rowsRic, p -> matDummySizeR_struct, 
					   p -> matDummySizeR_2_struct, p-> vecDummySizeR2_struct);
	m = 0;
	for ( k=0; k<rowsRic; k++)
		for (j=0; j<rowsRic; j++, m++)
			f[m+position_R_y_Riccati] = gsl_matrix_get(p -> matdR_struct,k,j);
	/********************************************************************************
	 ***    If I am integrating the eigenfunctions V and Rij are interpolated,
	 ***    then perform linear interpolation: R contains the interpolated values
	 ********************************************************************************/				
	if (integratingV && interpolate_rij)
	{
		lowEndCsiCalc = p->csi_calc_struct[p->startInterpol_rij_struct];
		highEndCsiCalc = p->csi_calc_struct[p->endInterpol_rij_struct];
		if (((integratingV_fitSurface && (csi < lowEndCsiCalc || csi > highEndCsiCalc)))
			|| ((integratingV_fitCenter && (csi > lowEndCsiCalc || csi < highEndCsiCalc))))
		{			
			errorMessageAndExit("gsl_odeiv_evolve_Riccati.c -- func_Riccati_R", 
								"csi outside allowed range for interpolation of r_ij");			
		}
		else
			calculateInterpolated_rij(p->startInterpol_rij_struct, p->endInterpol_rij_struct, csi, p->csi_calc_struct, 							
									  p->Rij_calc_struct, p->matR_struct, rowsRic);		
		

	}//if (integratingV && interpolate_rij)	
	/********************************************************************************
	 ***    If I am integrating R only, set dV/dx = dR/dx for each row of Riccati
	 ********************************************************************************/				
	for (i=0; i < rowsRic; i++)
		f[i + position_V_y_Riccati] = f[i + position_R_y_Riccati];	
	/********************************************************************************
	 ***   If I am integrating V then calculate the proper terms
	 ********************************************************************************/				
	if (integratingV)
	{
		/* perform CR + D and store it in the dummy matrix */
		matrix_mul_rows_by_cols(p->matC_struct, p->matR_struct, p->matCR_struct, rowsRic, 
								p->vecDummySizeR2_struct);
		matrix_sum (p->matCR_struct, p->matD_struct, p->matDummySizeR_struct, rowsRic);
		
		/* (CR + D)* vector */
		for (i=0; i < rowsRic; i++)
			gsl_vector_set(p->vecDummySizeR_struct, i, y[i + position_V_y_Riccati]);
		
		mul_matrix_by_vec(p->matDummySizeR_struct, p->vecDummySizeR_struct, p->vecDummySizeR_struct, p->vecDummySizeR_3_struct, rowsRic);
		
		/*dVi*/
		for (i=0; i < rowsRic; i++)
			f[i + position_V_y_Riccati] = gsl_vector_get(p->vecDummySizeR_struct,i);
		
		if(interpolate_rij)
		{
			for (i = 0; i <rowsRic; i++)
				for (j = 0; j <rowsRic; j++)
					f[rowsRic*i+j+position_R_y_Riccati] = f[j+position_V_y_Riccati];
		}
	}//if (integratingV)
		
	return GSL_SUCCESS;
}

int jac_Riccati_R (double csi, const double y[], double *dfdy, double dfdt[], void *params)
{
	/**************************************************************************
	 ************           
	 ************    Integrator for eigenfunctions and eigenfrequencies
	 ************    following Takata & Loffler 2004
	 ************              2004PASJ...56..645T
	 ************           
	 ************    See Eq.(13) and Eq. (22)
	 ************     
	 ************     dR/dx = B+AR-RD-RCR
	 ************     dv/dx = (CR+D)v
	 ************     
	 ************     The Jacobian terms are:
	 ************     
	 ************     d/dR(dR/dx) = A*dR/dR - dR/dR*D - dR/dR * CR - RC*dR/dR
	 ************     d/v(dR/dx) = 0
	 ************     d/dR(dv/dx) = C*dR/dR*v
	 ************     d/v(dv/dx) = CR+D
	 ************     
	 ************     General configuration of the variables for R = 2x2:
	 ************     y[0] = r00........y[1] = r01
	 ************     y[2] = r10........y[3] = r11
	 ************     y[4] = v0.........y[5] = v1 	 
	 ************     
	 *************************************************************************/
	params_integrator_adiabatic_struct *p = (params_integrator_adiabatic_struct *) params;
	
	int i = 0, j = 0, k = 0, n = 0, t = 0;
	int rowsRic	         = p->rowsRic_struct, 
	sizeR                = rowsRic*rowsRic, 
	sizeJacobian         = sizeR+rowsRic; 
	const int position_R_y_Riccati = 0, 
	position_V_y_Riccati = sizeR;

	int integratingRonly = 0, integratingV = 0;
	
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, sizeJacobian, sizeJacobian);
    gsl_matrix * m = &dfdy_mat.matrix; 
	
	if(p->inWard_outWard_struct>1){integratingRonly = 1;}
	if(p->inWard_outWard_struct<=1){integratingV = 1;}
	
	gsl_matrix_set_zero(m);

	/*****************************************************************
	 *****     d/dR(dR/dx) = A*dR/dR - dR/dR*D - dR/dR * CR - RC*dR/dR
	 ******************************************************************/
	/* (A- RC)*dR/dR term above*/
	for(n = 0; n < rowsRic; n++)
		for(k = 0; k<rowsRic; k++)
		{
			j = 0;
			for(i = 0; i <=  rowsRic * (rowsRic-1); i = i+rowsRic, j++)		
			{
				gsl_matrix_set (m, k+n*rowsRic, i+k, gsl_matrix_get(m, k+n*rowsRic, i+k) + gsl_matrix_get(p -> matA_struct,n,j));
				
				for(t = 0; t < rowsRic; t++)
					gsl_matrix_set (m, k+n*rowsRic, i+k, gsl_matrix_get(m, k+n*rowsRic, i+k) - y[t+n*rowsRic] * gsl_matrix_get(p -> matC_struct,t,j));				
			}
		}	
	
	/* - dR/dR(D - CR) term above*/
	for (n = 0; n<rowsRic; n++)			
		for(i=0; i < rowsRic; i++)
			for(j = 0; j<rowsRic; j++)
			{
				gsl_matrix_set (m, i+n*rowsRic, j+n*rowsRic, gsl_matrix_get(m, i+n*rowsRic, j+n*rowsRic) - gsl_matrix_get(p -> matD_struct,j,i));
				for(k = 0; k < rowsRic; k++)	
					gsl_matrix_set (m, i+n*rowsRic, j+n*rowsRic, gsl_matrix_get(m, i+n*rowsRic, j+n*rowsRic) - gsl_matrix_get(p -> matC_struct,j,k)*y[k*rowsRic+i]);
			}
	
	
	/* d/dR(dv/dx) = C*dR/dR*v*/
	for (n = 0; n<rowsRic; n++)
		for(j=0; j < rowsRic; j++)
			for(i = 0; i<rowsRic; i++)
				gsl_matrix_set (m, rowsRic*rowsRic+n, i+j*rowsRic, gsl_matrix_get(m, rowsRic*rowsRic+n, i+j*rowsRic) 
								+ gsl_matrix_get(p -> matC_struct,n,j)*y[position_V_y_Riccati+i]);
	
	
	/*d/v(dv/dx) = CR+D*/
	for(j = 0; j < rowsRic; j++)
		for(i = 0; i < rowsRic; i++)
		{
			gsl_matrix_set(m, sizeR+j, sizeR+i, gsl_matrix_get(m, sizeR+j, sizeR+i) + gsl_matrix_get(p -> matD_struct,j,i));
			for(k = 0; k<rowsRic; k++)
				gsl_matrix_set(m, sizeR+j, sizeR+i, gsl_matrix_get(m, sizeR+j, sizeR+i) + gsl_matrix_get(p -> matC_struct,j,k)*y[i+k*rowsRic]);
		}
	
	/*****************************************************************
	 *****  If I am interpolating R during the calculation of V, then:
	 *****   y[0] = v0........y[1] = v1  
	 *****   y[2] = v0........y[3] = v1  
	 *****   y[4] = v0........y[5] = v1  
	 *****     
	 *****  d/dv0(dv0/dx), d/dv1(dv0/dx), d/dv0(dv1/dx), d/dv1(dv1/dx)   	 
	 *****  are the elements used to calculate the full Jacobian and 
  	 *****  the interpolated Rij are used
	 ******************************************************************/
	
	if(integratingV && p->interpol_rij_struct)
	{
		gsl_matrix_set_zero(m);
		
		/*d/v(dv/dx) in the right bottom corner of the Jacobian */
		for(j = 0; j < rowsRic; j++)
			for(i = 0; i < rowsRic; i++)
			{
				gsl_matrix_set(m, sizeR+j, sizeR+i, gsl_matrix_get(m, sizeR+j, sizeR+i) + gsl_matrix_get(p -> matD_struct,j,i));
				for(k = 0; k<rowsRic; k++)
					gsl_matrix_set(m, sizeR+j, sizeR+i, gsl_matrix_get(m, sizeR+j, sizeR+i) + gsl_matrix_get(p -> matC_struct,j,k)*gsl_matrix_get(p->matR_struct, k, i));
			}
		
		
		/* part of the Jacobian for d/R(dv/dx) copied from the one calculated above */
		for(k=0; k<rowsRic; k++)
			for(i=0; i<rowsRic; i++)
				for(j=0; j<rowsRic; j++)
					gsl_matrix_set(m, k + position_V_y_Riccati, j + i*rowsRic + position_R_y_Riccati, 
								   gsl_matrix_get(m, k+position_V_y_Riccati, j+position_V_y_Riccati));
		
		
		/* part of the Jacobian for d/R(dR/dx) and d/v(dR/dx) copied from the above */

		for(i=0; i<rowsRic; i++)
			for(j=0; j<rowsRic; j++)
				for(k = 0; k<sizeJacobian; k++)
					gsl_matrix_set(m, j + i*rowsRic + position_R_y_Riccati, k , gsl_matrix_get(m, j+position_V_y_Riccati, k));

	}//if(integratingV && p->interpol_rij_struct)
	
	/*****************************************************************
	 *****  If I am integrating R only then:
	 *****   y[0] = r00........y[1] = r01  
	 *****   y[2] = r10........y[3] = r11  
	 *****   y[3] = r00........y[4] = r01  
	 ******************************************************************/
	if(integratingRonly)
	{		
		/* part of the Jacobian for d/v(dR/dx) copied from d/R(dR/dx) */
		for(i=0; i<sizeR; i++)
			for(j=sizeR; j<sizeJacobian; j++)
				gsl_matrix_set(m,i ,j, gsl_matrix_get(m,i,j-sizeR));
		
		/* part of the Jacobian for d/v(dv/dx) and d/R(dv/dx) copied from d/v(dR/dx) and d/R(dR/dx) */
		for(i=sizeR; i<sizeR+rowsRic; i++)
			for(j=0; j<sizeJacobian; j++)
				gsl_matrix_set(m,i,j, gsl_matrix_get(m,i-sizeR,j));
		
	}//if(integratingRonly)

	
	for (i=0; i<sizeJacobian; i++){dfdt[i] = 0.0;}
		
	return GSL_SUCCESS;
}

/********************************************************************************
 ********************************************************************************
 ******************    
 ******************				non Adiabatic case:
 ******************    
 ********************************************************************************
 ********************************************************************************/					
int func_Riccati_Ronly (double csi, const double y[], double f[], void *params)
{	
	
	int j = 0, k = 0, m = 0;
  	params_integrator_R_nonAd_struct *p = (params_integrator_R_nonAd_struct *) params;

	/********************************************************************************
	 ***    Interpolating with Steffen the stellar model's quantities
	 ********************************************************************************/						
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVgCoeffs_struct,csi, p->VgFuncs_struct);	
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitAstarCoeffs_struct,csi, p->AstarFuncs_struct);
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitUcoeffs_struct,csi, p->UFuncs_struct);
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC1coeffs_struct,csi, p->c1Funcs_struct); 

	/*Terms entering non adiabatic equations*/
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVCoeffs_struct,csi, p->VFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitDelADcoeffs_struct,csi, p->delADFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitDelCoeffs_struct,csi, p->delFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVtCoeffs_struct,csi, p->VtFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitKsCoeffs_struct,csi, p->ksFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC2coeffs_struct,csi, p->c2Funcs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC4coeffs_struct,csi, p->c4Funcs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitdlnLR_dlnrCoeffs_struct,csi, p->dlnLR_dlnrFuncs_struct); 

	/*these should be very small for a WD*/
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitEpsADcoeffs_struct,csi, p->epsADFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitEpsScoeffs_struct,csi, p->epsSFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC3coeffs_struct,csi, p->c3Funcs_struct); 
	/********************************************************************************
	 ***    Starting from original ABCD, going to M, permuting M and back to ABCD
	 ***    after this:T and M are unchanged, while A, B, C, D are permuted 
	 ********************************************************************************/				
	permute_Riccati_ABCD(csi, p->VgFuncs_struct[0], p->AstarFuncs_struct[0], 
						 p->UFuncs_struct[0], p->c1Funcs_struct[0], p->l_struct, p->omega_r_struct, 						 
						 p->omega_i_struct, p->VtFuncs_struct[0], p->VFuncs_struct[0], 
						 p->delFuncs_struct[0], p->delADFuncs_struct[0], p->ksFuncs_struct[0], 
						 p->c2Funcs_struct[0], p->c3Funcs_struct[0], p->c4Funcs_struct[0], 
						 p->epsADFuncs_struct[0], p->dlnLR_dlnrFuncs_struct[0], p->epsSFuncs_struct[0], 
						 p->matA_struct, p->matB_struct, p->matC_struct,p->matD_struct, 
						 p->matM_struct, p->matT_struct,
						 p-> matDummySizeT_struct, p-> matDummySizeT_2_struct, 
						 p->vecDummySizeT2_struct, p->permDummySizeT_struct, p->rowsRic_struct, 
						 p->inWard_outWard_struct, p->nonAdiabatic_struct, p->tidesFlag_struct);
	
	/* ABCD are simmetric by construction*/
	/********************************************************************************
	 ***  Fill in Riccati matrix to create Riccati equations ==>
	 ***  Matrix R contains the calculated values
	 ********************************************************************************/					
	m=0;
	for (k=0; k<static_cast<int>(p->rowsRic_struct); k++)
		for (j=0; j<p->rowsRic_struct; j++, m++)
			gsl_matrix_set(p->matR_struct,k,j, y[m]);	
	/********************************************************************************
	 ***    Create Riccati equations:
	 ***    after this:A, B, C, D and R are unchanged, while I have dR/dx
	 ********************************************************************************/				
	create_Riccati_eqs(p -> matA_struct, p -> matB_struct, p -> matC_struct, p -> matD_struct, 
					   p -> matR_struct, p -> matdR_struct, p->rowsRic_struct, p -> matDummySizeR_struct, 
					   p -> matDummySizeR_2_struct, p-> vecDummySizeR2_struct);
	m = 0;
	for ( k=0; k<p->rowsRic_struct; k++)
		for (j=0; j<p->rowsRic_struct; j++, m++)
			f[m] = gsl_matrix_get(p -> matdR_struct,k,j);
	
	return GSL_SUCCESS;
}


int jac_Riccati_Ronly (double csi, const double y[], double *dfdy, double dfdt[], void *params)
{
	/**************************************************************************
	 ************           
	 ************    Integrator for eigenfrequencies
	 ************    following Takata & Loffler 2004
	 ************              2004PASJ...56..645T
	 ************           
	 ************    See Eq.(13):
	 ************     
	 ************     dR/dx = B+AR-RD-RCR
	 ************     
	 ************     The Jacobian terms are:
	 ************     
	 ************     d/dR(dR/dx) = A*dR/dR - dR/dR*D - dR/dR * CR - RC*dR/dR
	 ************     
	 ************     General configuration of the variables for R = 6x6:
	 ************     y[0] = r00........y[5] = r05
	 ************     y[6] = r10........y[11] = r15
	 ************     y[12] = r20.......y[17] = r25
	 ************     y[18] = r30.......y[23] = r35
	 ************     y[24] = r40.......y[29] = r45
	 ************     y[30] = r50.......y[35] = r55
	 ************     
	 *************************************************************************/
	params_integrator_R_nonAd_struct *p = (params_integrator_R_nonAd_struct *) params;
	
	int i = 0, k = 0, n = 0, j = 0, t = 0,
	sizeJacobian = p->rowsRic_struct*p->rowsRic_struct; 

	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, sizeJacobian, sizeJacobian);
    gsl_matrix * m = &dfdy_mat.matrix; 
	
	gsl_matrix_set_zero(m);
	
	/* (A- RC)*dR/dR term above*/
	for(n = 0; n <  p->rowsRic_struct; n++)
		for(k = 0; k< p->rowsRic_struct; k++)
		{
			j = 0;
			for(i = 0; i <=   p->rowsRic_struct * ( p->rowsRic_struct-1); i = i+ p->rowsRic_struct, j++)		
			{
				gsl_matrix_set (m, k+n* p->rowsRic_struct, i+k, gsl_matrix_get(m, k+n* p->rowsRic_struct, i+k) + gsl_matrix_get(p -> matA_struct,n,j));
				
				for(t = 0; t <  p->rowsRic_struct; t++)
					gsl_matrix_set (m, k+n* p->rowsRic_struct, i+k, gsl_matrix_get(m, k+n* p->rowsRic_struct, i+k) - y[t+n* p->rowsRic_struct] * gsl_matrix_get(p -> matC_struct,t,j));				
			}
		}	
	
	/* - dR/dR(D + CR) term above*/
	for (n = 0; n< p->rowsRic_struct; n++)			
		for(i=0; i <  p->rowsRic_struct; i++)
			for(j = 0; j< p->rowsRic_struct; j++)
			{
				gsl_matrix_set (m, i+n* p->rowsRic_struct, j+n* p->rowsRic_struct, 
								gsl_matrix_get(m, i+n* p->rowsRic_struct, j+n* p->rowsRic_struct) - gsl_matrix_get(p -> matD_struct,j,i));
				for(k = 0; k <  p->rowsRic_struct; k++)	
					gsl_matrix_set (m, i+n* p->rowsRic_struct, j+n* p->rowsRic_struct, gsl_matrix_get(m, i+n* p->rowsRic_struct, j+n* p->rowsRic_struct) - 
									gsl_matrix_get(p -> matC_struct,j,k)*y[k* p->rowsRic_struct+i]);
			}
	
	
	for (i=0; i<sizeJacobian; i++){dfdt[i] = 0.0;}
	
	return GSL_SUCCESS;
}

/********************************************************************************
 ***    Using rk4 or bsimp and interpolating with Steffen
 ********************************************************************************/				
int func_Riccati_Vonly (double csi, const double y[], double f[], void *params)
{	
	int i = 0;
  	params_integrator_V_nonAd_struct *p = (params_integrator_V_nonAd_struct *) params;
	
	double lowEndCsiCalc = 0.0, highEndCsiCalc = 0.0;
	/********************************************************************************
	 ***    forcing csi to be within integration interval
	 ***    to avoid numerical errors to the 16th decimal place when interpolating Rij
	 ********************************************************************************/				
	/*integrating fit to surface*/	
	if (p->inWard_outWard_struct == 0)
	{
		if (csi > p->csiFin_force_struct){csi = p->csiFin_force_struct;}
		if (csi < p->csiIn_force_struct){csi = p->csiIn_force_struct;}
	}
	/*integrating fit to center*/
	if (p->inWard_outWard_struct == 1)
	{
		if (csi < p->csiFin_force_struct){csi = p->csiFin_force_struct;}
		if (csi > p->csiIn_force_struct){csi = p->csiIn_force_struct;}
	}
	/********************************************************************************
	 ***    Interpolating with Steffen the stellar model's quantities
	 ********************************************************************************/				
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVgCoeffs_struct,csi, p->VgFuncs_struct);	
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitAstarCoeffs_struct,csi, p->AstarFuncs_struct);
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitUcoeffs_struct,csi, p->UFuncs_struct);
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC1coeffs_struct,csi, p->c1Funcs_struct); 

	/*Terms entering non adiabatic equations*/
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVCoeffs_struct,csi, p->VFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitDelADcoeffs_struct,csi, p->delADFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitDelCoeffs_struct,csi, p->delFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVtCoeffs_struct,csi, p->VtFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitKsCoeffs_struct,csi, p->ksFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC2coeffs_struct,csi, p->c2Funcs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC4coeffs_struct,csi, p->c4Funcs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitdlnLR_dlnrCoeffs_struct,csi, p->dlnLR_dlnrFuncs_struct); 

	/*these should be very small for a WD*/
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitEpsADcoeffs_struct,csi, p->epsADFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitEpsScoeffs_struct,csi, p->epsSFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC3coeffs_struct,csi, p->c3Funcs_struct); 
	/********************************************************************************
	 ***    Starting from original ABCD, going to M, permuting M and back to ABCD
	 ***    after this:T and M are unchanged, while A, B, C, D are permuted 
	 ********************************************************************************/				
	permute_Riccati_ABCD(csi, p->VgFuncs_struct[0], p->AstarFuncs_struct[0], 
						 p->UFuncs_struct[0], p->c1Funcs_struct[0], p->l_struct, p->omega_r_struct, 						 
						 p->omega_i_struct, p->VtFuncs_struct[0], p->VFuncs_struct[0], 
						 p->delFuncs_struct[0], p->delADFuncs_struct[0], p->ksFuncs_struct[0], 
						 p->c2Funcs_struct[0], p->c3Funcs_struct[0], p->c4Funcs_struct[0], 
						 p->epsADFuncs_struct[0], p->dlnLR_dlnrFuncs_struct[0], p->epsSFuncs_struct[0], 
						 p->matA_struct, p->matB_struct, p->matC_struct,p->matD_struct, 
						 p->matM_struct, p->matT_struct,
						 p-> matDummySizeT_struct, p-> matDummySizeT_2_struct, 
						 p->vecDummySizeT2_struct, p->permDummySizeT_struct, p->rowsRic_struct, 
						 p->inWard_outWard_struct, p->nonAdiabatic_struct, p->tidesFlag_struct);
	/********************************************************************************
	 ***    Interpolating Rij and create matrix R
	 ********************************************************************************/				
	lowEndCsiCalc = p->csi_calc_struct[p->startInterpol_rij_struct];	
	highEndCsiCalc = p->csi_calc_struct[p->endInterpol_rij_struct];
	
	if ((((p->inWard_outWard_struct == 0) && (csi < lowEndCsiCalc || csi > highEndCsiCalc)))
		|| (((p->inWard_outWard_struct == 1) && (csi > lowEndCsiCalc || csi < highEndCsiCalc))))
		errorMessageAndExit("gsl_odeiv_evolve_Riccati.c", "csi outside allowed range for interpolation of r_ij");			
	else
		evalSteffenInterp_RiccatiElmnts(p->sizeRij_calc_struct, p->startInterpol_rij_struct, p->endInterpol_rij_struct, 
										p->csi_calc_struct, p->fitRij_calcCoeffs_struct, csi, p->Rij_calcFuncs_struct, 
										p->matR_struct, p->rowsRic_struct);
	/********************************************************************************
	 ***    Create equations for dV/dx
	 ********************************************************************************/				
	/* perform CR + D and store it in the dummy matrix */
	matrix_mul_rows_by_cols(p->matC_struct, p->matR_struct, p->matCR_struct, p->rowsRic_struct, 
							p->vecDummySizeR2_struct);
	matrix_sum (p->matCR_struct, p->matD_struct, p->matDummySizeR_struct, p->rowsRic_struct);
	
	/* (CR + D)* vector */
	for (i=0; i < p->rowsRic_struct; i++)
		gsl_vector_set(p->vecDummySizeR_struct, i, y[i]);
	
	mul_matrix_by_vec(p->matDummySizeR_struct, p->vecDummySizeR_struct, p->vecDummySizeR_struct, p->vecDummySizeR_3_struct, p->rowsRic_struct);
	
	/*dVi*/
	for (i=0; i < p->rowsRic_struct; i++)
		f[i] = gsl_vector_get(p->vecDummySizeR_struct,i);

	return GSL_SUCCESS;
}
/********************************************************************************
 ***    Using rkf45 and interpolating using the integrator state
 ********************************************************************************/				
int func_Riccati_Vonly_rkf45 (double csi, const double y[], double f[], void *params)
{	
	int i = 0;
  	params_integrator_V_nonAd_rkf45_struct *p = (params_integrator_V_nonAd_rkf45_struct *) params;
	
	double lowEndCsiCalc = 0.0, highEndCsiCalc = 0.0;
	/********************************************************************************
	 ***    forcing csi to be within integration interval
	 ***    to avoid numerical errors to the 16th decimal place when interpolating Rij
	 ********************************************************************************/				
	/*integrating fit to surface*/
	if (p->inWard_outWard_struct == 0)
	{
		if (csi > p->csiFin_force_struct){csi = p->csiFin_force_struct;}
		if (csi < p->csiIn_force_struct){csi = p->csiIn_force_struct;}
	}
	/*integrating fit to center*/
	if (p->inWard_outWard_struct == 1)
	{
		if (csi < p->csiFin_force_struct){csi = p->csiFin_force_struct;}
		if (csi > p->csiIn_force_struct){csi = p->csiIn_force_struct;}
	}
	/********************************************************************************
	 ***    Interpolating with Steffen the stellar model's quantities
	 ********************************************************************************/				
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVgCoeffs_struct,csi, p->VgFuncs_struct);	
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitAstarCoeffs_struct,csi, p->AstarFuncs_struct);
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitUcoeffs_struct,csi, p->UFuncs_struct);
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC1coeffs_struct,csi, p->c1Funcs_struct); 
	
	/*Terms entering non adiabatic equations*/
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVCoeffs_struct,csi, p->VFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitDelADcoeffs_struct,csi, p->delADFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitDelCoeffs_struct,csi, p->delFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitVtCoeffs_struct,csi, p->VtFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitKsCoeffs_struct,csi, p->ksFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC2coeffs_struct,csi, p->c2Funcs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC4coeffs_struct,csi, p->c4Funcs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitdlnLR_dlnrCoeffs_struct,csi, p->dlnLR_dlnrFuncs_struct); 
	
	/*these should be very small for a WD*/
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitEpsADcoeffs_struct,csi, p->epsADFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitEpsScoeffs_struct,csi, p->epsSFuncs_struct); 
	evalSteffenInterp(p->polyMesh_struct,p->csiRel_vec_struct,p->fitC3coeffs_struct,csi, p->c3Funcs_struct); 
	/********************************************************************************
	 ***    Starting from original ABCD, going to M, permuting M and back to ABCD
	 ***    after this:T and M are unchanged, while A, B, C, D are permuted 
	 ********************************************************************************/				
	permute_Riccati_ABCD(csi, p->VgFuncs_struct[0], p->AstarFuncs_struct[0], 
						 p->UFuncs_struct[0], p->c1Funcs_struct[0], p->l_struct, p->omega_r_struct, 						 
						 p->omega_i_struct, p->VtFuncs_struct[0], p->VFuncs_struct[0], 
						 p->delFuncs_struct[0], p->delADFuncs_struct[0], p->ksFuncs_struct[0], 
						 p->c2Funcs_struct[0], p->c3Funcs_struct[0], p->c4Funcs_struct[0], 
						 p->epsADFuncs_struct[0], p->dlnLR_dlnrFuncs_struct[0], p->epsSFuncs_struct[0], 
						 p->matA_struct, p->matB_struct, p->matC_struct,p->matD_struct, 
						 p->matM_struct, p->matT_struct,
						 p-> matDummySizeT_struct, p-> matDummySizeT_2_struct, 
						 p->vecDummySizeT2_struct, p->permDummySizeT_struct, p->rowsRic_struct, 
						 p->inWard_outWard_struct, p->nonAdiabatic_struct, p->tidesFlag_struct);
	/********************************************************************************
	 ***    Interpolating Rij and create matrix R
	 ********************************************************************************/				
	lowEndCsiCalc = p->csi_calc_struct[p->startInterpol_rij_struct];	
	highEndCsiCalc = p->csi_calc_struct[p->endInterpol_rij_struct];
	
	if ((((p->inWard_outWard_struct == 0) && (csi < lowEndCsiCalc || csi > highEndCsiCalc)))
		|| (((p->inWard_outWard_struct == 1) && (csi > lowEndCsiCalc || csi < highEndCsiCalc))))
		errorMessageAndExit("gsl_odeiv_evolve_Riccati.c", "csi outside allowed range for interpolation of r_ij");			
	else
			computing_Rij_from_rkf45_state(p->Rij_calc_struct, p->rRel_and_rkf45_state_struct, p->Rij_from_rkf45_state_struct,
										   p->startInterpol_rij_struct, p->endInterpol_rij_struct, csi, p->matR_struct, p->rowsRic_struct, 
										   p->sizeRij_calc_struct);
	
	/********************************************************************************
	 ***    Create equations for dV/dx
	 ********************************************************************************/				
	/* perform CR + D and store it in the dummy matrix */
	matrix_mul_rows_by_cols(p->matC_struct, p->matR_struct, p->matCR_struct, p->rowsRic_struct, 
							p->vecDummySizeR2_struct);
	matrix_sum (p->matCR_struct, p->matD_struct, p->matDummySizeR_struct, p->rowsRic_struct);
	
	/* (CR + D)* vector */
	for (i=0; i < p->rowsRic_struct; i++)
		gsl_vector_set(p->vecDummySizeR_struct, i, y[i]);
	
	mul_matrix_by_vec(p->matDummySizeR_struct, p->vecDummySizeR_struct, p->vecDummySizeR_struct, p->vecDummySizeR_3_struct, p->rowsRic_struct);
	
	/*dVi*/
	for (i=0; i < p->rowsRic_struct; i++)
		f[i] = gsl_vector_get(p->vecDummySizeR_struct,i);
	
	return GSL_SUCCESS;
}

int jac_Riccati_Vonly (double csi, const double y[], double *dfdy, double dfdt[], void *params)
{
	/**************************************************************************
	 ************           
	 ************    Integrator for eigenfunctions and eigenfrequencies
	 ************    following Takata & Loffler 2004
	 ************              2004PASJ...56..645T
	 ************           
	 ************    Eq. (22)
	 ************     
	 ************     dv/dx = (CR+D)v
	 ************     
	 ************     The Jacobian terms are:
	 ************     
	 ************     d/v(dv/dx) = CR+D
	 ************     
	 ************     General configuration of the variables for R = 6x6:
	 ************     y[0] = v0
	 ************     y[1] = v1
	 ************     y[2] = v2
	 ************     y[3] = v3
	 ************     y[4] = v4
	 ************     y[5] = v5
	 ************     
	 *************************************************************************/
	params_integrator_V_nonAd_struct *p = (params_integrator_V_nonAd_struct *) params;
	
	int i = 0, j = 0, 
	sizeJacobian         = p->rowsRic_struct; 
	
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, sizeJacobian, sizeJacobian);
    gsl_matrix * m = &dfdy_mat.matrix; 
		
	gsl_matrix_set_zero(m);
	
	/*d/v(dv/dx) = CR+D*/
	for(j = 0; j < p->rowsRic_struct; j++)
		for(i = 0; i < p->rowsRic_struct; i++)
			gsl_matrix_set(m, j, i, gsl_matrix_get(m, j, i) + gsl_matrix_get(p->matCR_struct, j, i) + gsl_matrix_get(p -> matD_struct,j,i));
	
	for (i=0; i<sizeJacobian; i++){dfdt[i] = 0.0;}
	
	return GSL_SUCCESS;
}


