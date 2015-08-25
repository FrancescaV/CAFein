typedef struct
{
	double &csiIn_force_struct, &csiFin_force_struct;
	int rowsRic_struct, polyMesh_struct;
	
	double &omega_r_struct, &omega_i_struct, l_struct;
	
	double *VgFuncs_struct, *AstarFuncs_struct, *UFuncs_struct, *c1Funcs_struct, 
	*VFuncs_struct,  *delADFuncs_struct,  *delFuncs_struct,  *VtFuncs_struct,  *ksFuncs_struct,  *epsADFuncs_struct,  
	*epsSFuncs_struct,  *c2Funcs_struct,  *c3Funcs_struct,  *c4Funcs_struct,  *dlnLR_dlnrFuncs_struct; 
	
	double **fitVgCoeffs_struct, **fitAstarCoeffs_struct, **fitUcoeffs_struct, **fitC1coeffs_struct, 
	*csiRel_vec_struct;
	
	gsl_matrix *matA_struct,  *matB_struct, *matC_struct, *matD_struct, *matCR_struct,
	*matR_struct, *matdR_struct,
	*matT_struct, *matM_struct, *matDummySizeR_struct, *matDummySizeR_2_struct, 
	*matDummySizeT_struct, *matDummySizeT_2_struct;
	
	gsl_vector *vecDummySizeR2_struct, *vecDummySizeT2_struct, *vecDummySizeR_struct, *vecDummySizeR_3_struct;
	
	gsl_permutation *permDummySizeT_struct;
	
	int &inWard_outWard_struct, interpol_rij_struct;
	
	double **Rij_calc_struct, *csi_calc_struct;
		
	int &startInterpol_rij_struct, &endInterpol_rij_struct, 
	nonAdiabatic_struct, tidesFlag_struct, WD_tides_flag_func;

}params_integrator_adiabatic_struct;


