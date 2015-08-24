#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "vector"
#include "SteffenInterp.h"
#include "IOfiles.h"
#include "dma.h"
#include "numericalIntegration.h"
#include "rkf45_state.h"

using namespace std;

/*************************************************************************************
 ************ 
 ************	The routines described below are used to do some calculations on the
 ************	eigenfunctions or towards the calculation of the eigenfunctions.
 ************ 
 **************************************************************************************/


/*************************************************************************************
 ****** This routine is used to calculate the normalization of the eigenfunctions
 ****** according to Eq. 5.17 in Bart's thesis (equivalent to Eq. 9 in Fuller-Lai 2010)
 ****** 
 ****** NOTE--only for adiabatic--
 ******  
 **************************************************************************************/

double calcInnerProduct(int nSteps, int mesh_func1, int mesh_func2, double l, 
						double *csiRel_func1, double *xi_r_func1, double *xi_h_func1, double *rho_func1, 
						double *csiRel_func2, double *xi_r_func2, double *xi_h_func2, double &error)
{
	double *x, *f;
	x = dfun(nSteps);
	f = dfun(nSteps);

	double 
	**fitRho_Coeffs1, **fitXi_r_Coeffs1, **fitXi_h_Coeffs1, 
	*rho_Funcs1, *xi_r_Funcs1, *xi_h_Funcs1,
	**fitXi_r_Coeffs2, **fitXi_h_Coeffs2, 
	*xi_r_Funcs2, *xi_h_Funcs2;
	
	double innerProduct = 0.0, 
	dbl_i = 0.0, xMin = 0.0, xMax = 0.0, 
	dbl_n = static_cast<double>(nSteps);
	
	/* since I can not extrapolate with Steffen I have to integrate on the innermost interval */
	
	xMin = csiRel_func1[0];
	if(csiRel_func2[0] > csiRel_func1[0]){xMin = csiRel_func2[0];}

	xMax = csiRel_func1[mesh_func1-1];
	if(csiRel_func2[mesh_func2-1] < csiRel_func1[mesh_func1-1]){xMax = csiRel_func2[mesh_func2-1];}

	int i = 0;
	
	fitRho_Coeffs1  = dfunc(4, mesh_func1);
	fitXi_r_Coeffs1 = dfunc(4, mesh_func1);
	fitXi_h_Coeffs1 = dfunc(4, mesh_func1);
	fitXi_r_Coeffs2 = dfunc(4, mesh_func2);
	fitXi_h_Coeffs2 = dfunc(4, mesh_func2);
	
	rho_Funcs1  = dfun(4);
	xi_r_Funcs1 = dfun(4);
	xi_h_Funcs1 = dfun(4);
	xi_r_Funcs2 = dfun(4);
	xi_h_Funcs2 = dfun(4);
	
	calcSteffenInterp(mesh_func1, csiRel_func1, xi_r_func1, fitXi_r_Coeffs1);
	calcSteffenInterp(mesh_func1, csiRel_func1, xi_h_func1, fitXi_h_Coeffs1);
	calcSteffenInterp(mesh_func1, csiRel_func1, rho_func1,  fitRho_Coeffs1);

	calcSteffenInterp(mesh_func2, csiRel_func2, xi_r_func2, fitXi_r_Coeffs2);
	calcSteffenInterp(mesh_func2, csiRel_func2, xi_h_func2, fitXi_h_Coeffs2);
		
	for (i = 0; i<nSteps; i++)
	{
		dbl_i = static_cast<double>(i);
		x[i] = xMin + (dbl_i/dbl_n)*(xMax - xMin);
		
		evalSteffenInterp(mesh_func1,csiRel_func1,fitXi_r_Coeffs1,x[i],xi_r_Funcs1);
		evalSteffenInterp(mesh_func1,csiRel_func1,fitXi_h_Coeffs1,x[i],xi_h_Funcs1);
		evalSteffenInterp(mesh_func1,csiRel_func1,fitRho_Coeffs1,x[i],rho_Funcs1);

		evalSteffenInterp(mesh_func2,csiRel_func2,fitXi_r_Coeffs2,x[i],xi_r_Funcs2);
		evalSteffenInterp(mesh_func2,csiRel_func2,fitXi_h_Coeffs2,x[i],xi_h_Funcs2);
		
		f[i] = (xi_r_Funcs1[0]*xi_r_Funcs2[0]+l*(l+1.0)*(xi_h_Funcs1[0]*xi_h_Funcs2[0]))*pow(x[i],2.0)*rho_Funcs1[0];
	}
	innerProduct = FourPt_DoubleStar(x, f, nSteps, 0, nSteps-1, error);	

	
	free(x);
	free(f);
	
	free(rho_Funcs1);
	free(xi_r_Funcs1);
	free(xi_h_Funcs1);
	free(xi_r_Funcs2);
	free(xi_h_Funcs2);
	
	for (i=0; i<4; i++)
	{
		free(fitRho_Coeffs1[i]);
		free(fitXi_r_Coeffs1[i]);
		free(fitXi_h_Coeffs1[i]);
		free(fitXi_r_Coeffs2[i]);
		free(fitXi_h_Coeffs2[i]);
	}
	free(fitRho_Coeffs1);
	free(fitXi_r_Coeffs1);
	free(fitXi_h_Coeffs1);
	free(fitXi_r_Coeffs2);
	free(fitXi_h_Coeffs2);
	

	return innerProduct;
}

/*************************************************************************************
 ****** This routine is used to order the radial and orthogonal components of the
 ****** displacement (CsiR and CsiH, respectively). CsiR and CsiH are taken from 
 ****** the file eigenfunctions.dat, where they are calculated going from 
 ****** the fitting point to the star's surface and from the fitting point to
 ****** the star's center.
 ****** 
 ****** Input: 
 ****** - meshTot_eigenf_func = Total number of mesh points in the eigenfunction(s)
 ****** - elements_eigenf_func = matrix storing the calculated eigenfunctions
 ****** together with the radius
 ****** - rowsRic_func = # rows in the riccati matrix == # of eigenfunctions
 ****** 
 ****** Output:
 ****** - radiusEigenf_func = ordered radius for the eigenfunction(s)
 ****** - csi_Tr_func = vector containing the ordered radial component of the displacement(real and imaginary)
 ****** - csi_Th_func = vector containing the ordered orthogonal component of the displacement(real and imaginary)
 ****** 
 **************************************************************************************/

void orderRadius_CsiR_CsiH(int meshTot_eigenf_func, int meshTot_global_func, double *radiusEigenf_func,
						   double *csi_Tr_real_func, double *csi_Tr_imag_func, double *csi_Th_real_func, double *csi_Th_imag_func, 
						   vector<vector<double> > elements_eigenf_func, vector<vector<double> > elements_global_func, 
						   int rowsRic_func, int nonAdiabatic_func)
{
	
	int mesh_fitToS = 0, mesh_fitToC = 0, i = 0, j = 0; 
	double g = 0.0, y1r = 0.0, y2r = 0.0, y2i = 0.0, r = 0.0, w_r = 0.0, w_i = 0.0;

	/*real and imaginary omega*/
	w_r = elements_global_func[1][10];
	w_i = elements_global_func[1][11];
	/***************************************************************
	 ***    Counting number of mesh points for inward and outward 
	 ***    integration. Position of radius in the file doesn't 
	 ***    change between adiabatic and nonadiabatic configuration.
	 ***************************************************************/				
	for (i = 0; i< meshTot_eigenf_func; i++)
	{
		mesh_fitToS++;
		if (elements_eigenf_func[i+1][0] < elements_eigenf_func[i][0]){break;}
	}
	mesh_fitToC = meshTot_eigenf_func - mesh_fitToS;
	/***************************************************************
	 ***    Filling-in radius, csi_Tr, and csi_Th in an ordered way
	 ***    Note, the factor 6 below account for other columns that
	 ***	won't change if the size of R changes
	 ***************************************************************/				
	j = 0;
	for (i = meshTot_eigenf_func; i>=1; i--)
	{
		//old one
//		g = elements_eigenf_func[i][2+rowsRic_func+rowsRic_func+rowsRic_func*rowsRic_func];
		g = elements_eigenf_func[i][4+rowsRic_func*rowsRic_func];

		r = elements_eigenf_func[i][0];
		radiusEigenf_func[j] = r; 
		y1r = elements_eigenf_func[i][2];
		y2r = elements_eigenf_func[i][3];
		
		
//		if(nonAdiabatic_func)
//		{
//			y1i = elements_eigenf_func[i][8];
//			y2i = elements_eigenf_func[i][9];
//		}
//				
		csi_Tr_real_func[j] = y1r*r;
		csi_Th_real_func[j] = (g/((w_r*w_r + w_i*w_i)*(w_r*w_r + w_i*w_i)))*(y2r*(w_r*w_r - w_i*w_i) + 2.0*y2i*w_r*w_i);

//		if(nonAdiabatic_func)
//		{
//			csi_Tr_imag_func[j] = y1i*r;
//			csi_Th_imag_func[j] = (g/((w_r*w_r + w_i*w_i)*(w_r*w_r + w_i*w_i)))*(y2i*(w_r*w_r - w_i*w_i) - 2.0*y2r*w_r*w_i);
//		}
//		else
//		{
			csi_Tr_imag_func[j] = 0.0;
			csi_Th_imag_func[j] = 0.0;
//		}

		j++;
		if (j > mesh_fitToC){break;}

	
	}//for (i = meshTot_eigenf_func; i>=1; i--)

	for (i = 1; i<mesh_fitToS; i++)
	{
		g = elements_eigenf_func[i][4+rowsRic_func*rowsRic_func];
		
		r = elements_eigenf_func[i][0];
		radiusEigenf_func[j] = r; 		
		y1r = elements_eigenf_func[i][2];
		y2r = elements_eigenf_func[i][3];
		
//		if(nonAdiabatic_func)
//		{
//			y1i = elements_eigenf_func[i][8];
//			y2i = elements_eigenf_func[i][9];
//		}
		
		csi_Tr_real_func[j] = y1r*r;	
		csi_Th_real_func[j] = (g/((w_r*w_r + w_i*w_i)*(w_r*w_r + w_i*w_i)))*(y2r*(w_r*w_r - w_i*w_i) + 2.0*y2i*w_r*w_i);
		
//		if(nonAdiabatic_func)
//		{
//			csi_Tr_imag_func[j] = y1i*r;
//			csi_Th_imag_func[j] = (g/((w_r*w_r + w_i*w_i)*(w_r*w_r + w_i*w_i)))*(y2i*(w_r*w_r - w_i*w_i) - 2.0*y2r*w_r*w_i);
//		}
//		else
//		{
			csi_Tr_imag_func[j] = 0.0;
			csi_Th_imag_func[j] = 0.0;
//		}
		
		j++;
	}//for (i = 1; i<mesh_fitToS; i++)
}
/*************************************************************************************
 ****** This routine is used during the calculation of the eigenfunctions to calculate
 ****** the Riccati coefficients Rij via linear interpolation
 ****** 
 ****** 
 ****** Input: 
 ****** - startLoop_func = left index of the interval to consider during the interpolation
 ****** - endLoop_func = right index of the interval to consider during the interpolation
 ****** - x_func = x-position where Rij have to be calculated
 ****** - x_vector_func = x-component of the calculated Rij
 ****** - Rii_calc_func = matrix containing the calculated Rij
 ****** - rowsRic_func = # of rows in the Riccati matrix R
 ****** 
 ****** Output
 ****** - matR_func = matrix containing the interpolated Rij at position x_func.
 ******  
 **************************************************************************************/
void calculateInterpolated_rij(int startLoop_func, int endLoop_func, double x_func, double *x_vector_func, 
							   double **Rii_calc_func, gsl_matrix *matR_func, int rowsRic_func)
{
	double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, A = 0.0, B = 0.0;
	int k = 0, x_low_index = 0, x_high_index = 0, i = 0, j = 0;
	
	/***************************************************************
	 ***    find interval where to interpolate
	 ***************************************************************/				
	
	for (k = startLoop_func + 1; k <= endLoop_func; k++)
		if ((x_func >= x_vector_func[k-1] && x_func <= x_vector_func[k]) || 
			(x_func <= x_vector_func[k-1] && x_func >= x_vector_func[k]))
		{
			x_low_index = k-1;
			x_high_index = k;
			break;
		}
	
	x1 = x_vector_func[x_low_index];
	x2 = x_vector_func[x_high_index];
		
	/***************************************************************
	 ***    interpolate
	 ***************************************************************/				
	
	k = 0;
	for (i = 0; i < rowsRic_func; i++)
		for (j = 0; j < rowsRic_func; j++)
		{
			y1 = Rii_calc_func[x_low_index][k];
			y2 = Rii_calc_func[x_high_index][k];
		
			A = (y1-y2)/(x1-x2);
			B = y1 - A*x1;
			gsl_matrix_set(matR_func,i,j, A*x_func + B);		
						
			k++;
		}	
}
/*************************************************************************************
 ****** 
 ****** Same as 'calculateInterpolated_rij' but specific for the non-adiabatic case
 ****** in which R = 6x6 and it is composed of 4 matrices of size 3x3 (R00, R01, R10, R11)
 ****** such that:
 ****** 
 ****** R00 = R11 and R10 = -R01
 ****** 
 ****** The only difference with the routine above is the way in which R is filled
 **************************************************************************************/

void calculateInterpolated_rij_R6by6(int startLoop_func, int endLoop_func, double x_func, double *x_vector_func, 
							   double **Rii_calc_func, gsl_matrix *matR_func, int rowsRic_func)
{
	double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, A = 0.0, B = 0.0;
	int k = 0, x_low_index = 0, x_high_index = 0, i = 0, j = 0, halfR = rowsRic_func/2;
	
	for (k = startLoop_func + 1; k <= endLoop_func; k++)
		if ((x_func >= x_vector_func[k-1] && x_func <= x_vector_func[k]) || 
			(x_func <= x_vector_func[k-1] && x_func >= x_vector_func[k]))
		{
			x_low_index = k-1;
			x_high_index = k;
			break;
		}
	
	x1 = x_vector_func[x_low_index];
	x2 = x_vector_func[x_high_index];
	
	k = 0;
	for (i = 0; i < halfR; i++)
		for (j = 0; j < rowsRic_func; j++)
		{
			y1 = Rii_calc_func[x_low_index][k];
			y2 = Rii_calc_func[x_high_index][k];

			A = (y1-y2)/(x1-x2);
			B = y1 - A*x1;
			gsl_matrix_set(matR_func,i,j, A*x_func + B);		

			if(j < halfR)
				gsl_matrix_set(matR_func,i+halfR,j+halfR, A*x_func + B);
			else
				gsl_matrix_set(matR_func,i+halfR,j-halfR, -(A*x_func + B));
			k++;
		}		
}

/*************************************************************************************
 ****** 
 ****** Same as 'calculateInterpolated_rij_R6by6' but the bisection method
 ****** is used to find the interval where to interpolate.
 ****** 
 ****** The number of bisections to perform is high (10000), but the bisection
 ****** procedure is stopped when I am left with 5 elements to consider.
 ****** 
 **************************************************************************************/

void calculateInterpolated_rij_R6by6_bisection(int startLoop_func, int endLoop_func, double x_func, double *x_vector_func, 
									 double **Rii_calc_func, gsl_matrix *matR_func, int rowsRic_func)
{
	double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, A = 0.0, B = 0.0;
	int k = 0, x_low_index = 0, x_high_index = 0, i = 0, j = 0, difference = 0, split = 0, increasingX = 0,
	halfR = rowsRic_func/2, numberBisections = 10000;
		
	/***************************************************************
	 ***    Determining if csi is increasing or decreasing
	 ***   (i.e. if I am in the interval fit-surface or fit-center)	 
	 ***************************************************************/					
	if (x_vector_func[endLoop_func] > x_vector_func[startLoop_func]){increasingX = 1;}

	for (k=0; k<numberBisections; k++)
	{
		/***************************************************************
		 ***   Determining the distance between the beginning and the
		 ***   end of the interval	 
		 ***************************************************************/							
		difference = endLoop_func-startLoop_func;
		if (difference < 5) {break;}					

		if (difference%2 == 0) {;}
		else {difference = difference +1;}

		/***************************************************************
		 ***   splitting the interval in 2:
		 ***   
		 ***   If the searched root is on the left interval, I move the 
		 ***   right end of the right interval to the middle
		 ***   If the searched root is on the right interval, I move the 
		 ***   left end of the left interval to the middle.
		 ***************************************************************/					
		split = difference/2;

		if(increasingX)
		{
			if (x_func >= x_vector_func[startLoop_func] && x_func < x_vector_func[startLoop_func+split]) 
			{
				endLoop_func = startLoop_func+split; 
			}
			else {startLoop_func = startLoop_func+split;}
		}
		else
		{
			if (x_func <= x_vector_func[startLoop_func] && x_func > x_vector_func[startLoop_func+split]) 
			{
				endLoop_func = startLoop_func+split; 
			}
			else {startLoop_func = startLoop_func+split;}
		}		
	}//for (k=0; k<numberBisections; k++)
	
	/***************************************************************
	 ***   Determining the interval where to perform interpolation	 
	 ***************************************************************/								
	for (k = startLoop_func + 1; k <= endLoop_func; k++)
		if ((x_func >= x_vector_func[k-1] && x_func <= x_vector_func[k]) || 
			(x_func <= x_vector_func[k-1] && x_func >= x_vector_func[k]))
		{
			x_low_index = k-1;
			x_high_index = k;
			break;
		}
	
	x1 = x_vector_func[x_low_index];
	x2 = x_vector_func[x_high_index];
	
	/***************************************************************
	 ***   Interpolate and fill in R.
	 ***************************************************************/								
	k = 0;
	for (i = 0; i < halfR; i++)
		for (j = 0; j < rowsRic_func; j++)
		{
			y1 = Rii_calc_func[x_low_index][k];
			y2 = Rii_calc_func[x_high_index][k];
			
			A = (y1-y2)/(x1-x2);
			B = y1 - A*x1;
			gsl_matrix_set(matR_func,i,j, A*x_func + B);		
			
			if(j < halfR)
				gsl_matrix_set(matR_func,i+halfR,j+halfR, A*x_func + B);
			else
				gsl_matrix_set(matR_func,i+halfR,j-halfR, -(A*x_func + B));
			k++;
		}	
	
}
/*************************************************************************************
 ************ 
 ************	This routine takes as input file eigenfunctions.dat
 ************	and normalizes the eigenfunctions so that y8 = 1 at
 ************	the star's surface. Specifically: y8_real = 1 and y8_imag = 0
 ************ 
 ************	note: the number 2 refers to the starting position of yi in eigenfunctions.dat
 ************	
 ************	If y8 = y8Re + i*y8Im, I have to apply a rotation so that y8_real = 1 and y8_imag = 0
 ************	at r = R. This can be achieved if I multiply at the surface y8 by its complex conjugate and
 ************	then normalize it properly to make it of unitary lenght.
 ************	
 ************	So: y8_norm = |y8Norm|[cos(theta_S) + i sin(theta_S)]*[[cos(theta_S) - i sin(theta_S)]]
 ************	
 ************	
 ************	Given each eigenfunction y = yRe + i*yIm, the normalized eigenfunctions are given by:
 ************	
 ************	yReNorm = |yNorm|(cos(theta)cos(theta_S)+sin(theta)sin(theta_S))
 ************	yImNorm = |yNorm|(sin(theta)cos(theta_S)-sin(theta_S)cos(theta))
 ************	
 ************	where:
 ************	|yNorm| = |y|/|y8(r = R)|
 ************	theta = arg(y8(r = R))
 ************	
 **************************************************************************************/

void normalizeEigenfunctions_y8to1atR(const char *nameEigenfunctionsInput, const char *nameEigenfunctionsOutput, 
									  int adiabatic_func, int rowsRic_func, int tidesFlag_func, 
									  double arg_y8_atR_func, double modulus_y8_atR_func, int columnsEigenfuncFile_func)
{
	
	
	gsl_complex y_complex;
	/*file containing unnormalized eigenfunctions*/
	ifstream eigenfunctionsInput;
	eigenfunctionsInput.open(nameEigenfunctionsInput);
	
	/*file containing normalized eigenfunctions*/
 	ofstream eigenfunctionsOutput;
	eigenfunctionsOutput.open(nameEigenfunctionsOutput);
	eigenfunctionsOutput << setiosflags(ios_base::scientific)<< setprecision(16);

	//	labelF_eigenfunctions(eigenfunctionsOutput, rowsRic_func, adiabatic_func, tidesFlag_func);
	labelF_eigenfunctions_and_unperturbed_model(eigenfunctionsOutput, rowsRic_func, adiabatic_func, tidesFlag_func);
	
	int totalCol = 0, i = 0, elmntsOfUV = 0, numberOfRows = 0;	
	gsl_vector *eachRowInFile = gsl_vector_alloc(columnsEigenfuncFile_func+1);	
	
	/*adiabatic case: y1, y2, y3, y4, y7, y8 no imaginary part*/
	if(adiabatic_func){elmntsOfUV = rowsRic_func*2;}
	
	/*non adiabatic case: y1, y2, y3, y4, y5, y6, y7, y8 plus imaginary part*/
	else{elmntsOfUV = rowsRic_func;}
	
	double data = 0.0, yRe = 0.0, yIm = 0.0, yReNorm = 0.0, yImNorm = 0.0, 
	arg_y = 0.0, modulus_y = 0.0;

	while(eigenfunctionsInput)
	{
		string eachLine;		
		getline(eigenfunctionsInput, eachLine);
		stringstream elementsOfLine(eachLine);
		eachLine.clear();
		
		/*loop on each line filling in the vector eachRowInFile*/
		totalCol = 0;
		while(elementsOfLine)
		{
			elementsOfLine >> data;			
			gsl_vector_set(eachRowInFile, totalCol, data);			
			totalCol++;
		}	
		elementsOfLine.clear();
		
		/*Loop on each vector to normalize the eigenfunctions (skipping labels)*/
		if(totalCol > 1)
		{
			for(i = 2; i < elmntsOfUV+2; i++)
			{
				yRe = gsl_vector_get(eachRowInFile, i);
				
				if (adiabatic_func){yIm = 0.0;}

				else{yIm = gsl_vector_get(eachRowInFile, i+elmntsOfUV);}

				GSL_SET_COMPLEX(&y_complex, yRe, yIm);
				arg_y     = gsl_complex_arg(y_complex);
				modulus_y = gsl_complex_abs(y_complex);

				
				yReNorm = (modulus_y/modulus_y8_atR_func)*
				(cos(arg_y)*cos(arg_y8_atR_func)+sin(arg_y)*sin(arg_y8_atR_func));
				
				yImNorm = (modulus_y/modulus_y8_atR_func)*
				(sin(arg_y)*cos(arg_y8_atR_func)-sin(arg_y8_atR_func)*cos(arg_y));

				/*sobstitute the normalized eigenfunctions in the vector*/
				gsl_vector_set(eachRowInFile, i, yReNorm);
				if(adiabatic_func == 0){gsl_vector_set(eachRowInFile, i+elmntsOfUV, yImNorm);}
				
			}//for(i = 2; i < elmntsOfUV; i++)

			
			for(i = 0; i < totalCol-1; i++){eigenfunctionsOutput << gsl_vector_get(eachRowInFile, i)<< "  ";}				
			eigenfunctionsOutput << endl;
		}		
		numberOfRows++;
	}//while(eigenfunctionsInput)
	
	eigenfunctionsInput.close();
	eigenfunctionsOutput.close();
	
	gsl_vector_free(eachRowInFile);
}

/*************************************************************************************
 ****** This routine is used to order the radial components of the
 ****** displacement (CsiR) both for the dynamical tide and the static tide. 
 ****** Both CsiR_st and CsiR_dyn are taken from 
 ****** the file eigenfunctions.dat, 
 ****** 
 ****** Note:
 ******   
 ****** - csiR_dyn = y1*r
 ****** - csiR_st = -y3*r
 ******   
 ****** In eigenfunction.dat the variables are ordered going from 
 ****** the fitting point to the star's surface and from the fitting point to
 ****** the star's center.
 ****** 
 ****** Input: 
 ****** - meshTot_eigenf_func = Total number of mesh points in the eigenfunction(s)
 ****** - elements_eigenf_func = matrix storing the calculated eigenfunctions
 ****** together with the radius
 ****** 
 ****** Output:
 ****** - radiusEigenf_func = ordered radius for the eigenfunction(s)
 ****** - csi_Tr_dyn_func = vector containing the ordered radial component of the dynamic displacement(real and imaginary)
 ****** - csi_Tr_st_func = vector containing the ordered radial component of the static displacement(real and imaginary)
 ****** 
 **************************************************************************************/

void orderRadius_CsiR_linearityWithTides(int meshTot_eigenf_func, double *radiusEigenf_func, 
										 double *csi_Tr_real_dyn_func, double *csi_Tr_imag_dyn_func, 						   
										 double *csi_Tr_real_st_func, double *csi_Tr_imag_st_func, 						   
										 vector<vector<double> > elements_eigenf_func)
{
	
	int mesh_fitToS = 0, mesh_fitToC = 0, i = 0, j = 0; 
	double y1r = 0.0, y1i = 0.0, y3r = 0.0, y3i = 0.0, r = 0.0;
	
	/***************************************************************
	 ***    Counting number of mesh points for inward and outward 
	 ***    integration. Position of radius in the file doesn't 
	 ***    change between adiabatic and nonadiabatic configuration.
	 ***************************************************************/				
	for (i = 0; i< meshTot_eigenf_func; i++)
	{
		mesh_fitToS++;
		if (elements_eigenf_func[i+1][0] < elements_eigenf_func[i][0]){break;}
	}
	mesh_fitToC = meshTot_eigenf_func - mesh_fitToS;
	/***************************************************************
	 ***    Filling-in radius, csi_Tr, and csi_Th in an ordered way
	 ***    Note, the factor 6 below account for other columns that
	 ***	won't change if the size of R changes
	 ***************************************************************/				
	j = 0;
	for (i = meshTot_eigenf_func; i>=1; i--)
	{
		
		r = elements_eigenf_func[i][0];
		radiusEigenf_func[j] = r; 
		y1r = elements_eigenf_func[i][2];
		y1i = elements_eigenf_func[i][10];

		y3r = elements_eigenf_func[i][4];
		y3i = elements_eigenf_func[i][12];

		
		csi_Tr_real_dyn_func[j] = y1r*r;
		csi_Tr_imag_dyn_func[j] = y1i*r;

		csi_Tr_real_st_func[j] = -y3r*r;
		csi_Tr_imag_st_func[j] = -y3i*r;

		j++;
		if (j > mesh_fitToC){break;}				
	}//for (i = meshTot_eigenf_func; i>=1; i--)
	
	for (i = 1; i<mesh_fitToS; i++)
	{
		r = elements_eigenf_func[i][0];
		radiusEigenf_func[j] = r; 
		
		y1r = elements_eigenf_func[i][2];
		y1i = elements_eigenf_func[i][10];
		
		y3r = elements_eigenf_func[i][4];
		y3i = elements_eigenf_func[i][12];
				
		csi_Tr_real_dyn_func[j] = y1r*r;			
		csi_Tr_imag_dyn_func[j] = y1i*r;
		
		csi_Tr_real_st_func[j] = -y3r*r;
		csi_Tr_imag_st_func[j] = -y3i*r;
		
		j++;
	}//for (i = 1; i<mesh_fitToS; i++)
}



/*************************************************************************************
 ****** This routine is used to order the matrix containint eigenfunctions.dat, 
 ****** where the elements are calculated going from 
 ****** the fitting point to the star's surface and from the fitting point to
 ****** the star's center.
 **************************************************************************************/

void orderEigenfunctionElmnts(int meshTot_eigenf_func, vector<vector<double> > elements_eigenf_func, 
							  vector<vector<double> > &elements_eigenf_ordered_func)
{
	
	
	int mesh_fitToS = 0, mesh_fitToC = 0, i = 0, j = 0, k = 0; 
	double data = 0.0;
	
	vector<double> eachRow;

	int sizeCols = elements_eigenf_func[meshTot_eigenf_func].size();

	/***************************************************************
	 ***    Counting number of mesh points for inward and outward 
	 ***    integration. 
	 ***************************************************************/				
	for (i = 0; i< meshTot_eigenf_func; i++)
	{
		mesh_fitToS++;
		if (elements_eigenf_func[i+1][0] < elements_eigenf_func[i][0]){break;}
	}
	mesh_fitToC = meshTot_eigenf_func - mesh_fitToS;
	/***************************************************************
	 ***    Filling-in the new matrix with elements from center to
	 ***    surface.
	 ***************************************************************/				
	j = 0;

	for (i = meshTot_eigenf_func; i>=1; i--)
	{
		for (k = 0; k<sizeCols-1; k++)
		{
			data = elements_eigenf_func[i][k];
			eachRow.push_back(data);
		}		
		elements_eigenf_ordered_func.push_back(eachRow);
		eachRow.clear();	
		j++;
		if (j > mesh_fitToC){break;}
	}//for (i = meshTot_eigenf_func; i>=1; i--)
	
	for (i = 1; i<mesh_fitToS; i++)
	{
		for (k = 0; k<sizeCols-1; k++)
		{
			data = elements_eigenf_func[i][k];
			eachRow.push_back(data);
		}		
		elements_eigenf_ordered_func.push_back(eachRow);
		eachRow.clear();	
	}//for (i = 1; i<mesh_fitToS; i++)

}

/*************************************************************************************
 ****** The two routines below are used to store rkf45 state.
 ****** From rkf45's state, we can compute the elements Rij during the calculation of the
 ****** eigenfunctions.
 **************************************************************************************/
void storeIntegratorState_out(gsl_vector *vecIntegratorState_func, gsl_odeiv_step *step, double csiRel_func, int sizeRij_calc_func)
{
	
	rkf45_state_t *s_Riccati_state;
	s_Riccati_state = (rkf45_state_t *) step->state;
	
	gsl_vector_set(vecIntegratorState_func, 0, csiRel_func);
	
	int i = 0, j = 0;
	
	j = 1;
	for(i = 0; i < sizeRij_calc_func; i++)
	{
		gsl_vector_set(vecIntegratorState_func, j, s_Riccati_state->k1[i]);
		j++;
		gsl_vector_set(vecIntegratorState_func, j, s_Riccati_state->k6[i]);
		j++;
	}
	
}

void storeIntegratorState_in(gsl_matrix *matIntegratorState_func, gsl_odeiv_step *step, double csiRel_func, int sizeRij_calc_func, int position)
{
	
	rkf45_state_t *s_Riccati_state;
	s_Riccati_state = (rkf45_state_t *) step->state;
	
	gsl_matrix_set(matIntegratorState_func, position, 0, csiRel_func);
	
	int i = 0, j = 0;
	
	j = 1;
	for(i = 0; i < sizeRij_calc_func; i++)
	{
		gsl_matrix_set(matIntegratorState_func, position, j, s_Riccati_state->k1[i]);
		j++;
		gsl_matrix_set(matIntegratorState_func, position, j, s_Riccati_state->k6[i]);
		j++;
	}
	
}

/*************************************************************************************
 ****** Computing Rij using rkf45 state. In this routine I want to
 ****** use bysection to find the interval across which I interpolate Rij
 ******
 ****** Note: the zeroth column of the 2D vector rRel_and_rkf45_state_func contains
 ****** the star's radial profile
 **************************************************************************************/

void computing_Rij_from_rkf45_state(double **Rij_calc_func, double **rRel_and_rkf45_state_func, 
									double *Rij_from_rkf45_state, int position_xIn, int position_xFin, 
									double rRelStep, gsl_matrix *matR_func, int rowsRic_func, 
									int sizeRij_calc_func)
{
	int i = 0, k = 0, j = 0, r = 0, rOld = 0, increasingR = 0;
	double hUsed = 0.0, theta = 0.0, i0 = 0.0, i1 = 0.0, i6 = 0.0, iend = 0.0;
	
	/*checking if the radius is increasing or decreasing*/
	if (rRel_and_rkf45_state_func[position_xFin][0] > rRel_and_rkf45_state_func[position_xIn][0])
		increasingR = 1;
	
	/* Check if x0 is in the range of the (x_i) data set*/
	
	if (increasingR)
	{
		if (rRelStep < rRel_and_rkf45_state_func[position_xIn][0] || rRelStep > rRel_and_rkf45_state_func[position_xFin][0])
		{
			cout << "Extrapolation not allowed in " << "computing_Rij_from_rkf45_state!'" << endl;	         	
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		if (rRelStep > rRel_and_rkf45_state_func[position_xIn][0] || rRelStep < rRel_and_rkf45_state_func[position_xFin][0])
		{
			cout << "Extrapolation not allowed in " << "computing_Rij_from_rkf45_state!'" << endl;	         	
			exit(EXIT_FAILURE);
		}
	}
	
	/* Locate rRel in the rRel_and_rkf45_state_func[i][0] data set via bysection*/
	int xHighIndex = position_xFin, xLowIndex = position_xIn, difference = 0, split = 0;	
	int numberBisections = 10000;
	
	for (k=0; k<numberBisections; k++)
	{
		/* determining the interval where the crossing happens by bisection*/
		difference = xHighIndex-xLowIndex;
		if (difference < 10) {break;}					
		
		if (difference%2 == 0) {;}
		else {difference = difference +1;}
		
		split = difference/2;
		
		if(increasingR)
		{
			if (rRelStep >= rRel_and_rkf45_state_func[xLowIndex][0] && rRelStep < rRel_and_rkf45_state_func[xLowIndex+split][0]) 			
				xHighIndex = xLowIndex+split; 
			else {xLowIndex = xLowIndex+split;}
		}			
		else
		{
			if (rRelStep <= rRel_and_rkf45_state_func[xLowIndex][0] && rRelStep > rRel_and_rkf45_state_func[xLowIndex+split][0])
				xHighIndex = xLowIndex+split; 
			else {xLowIndex = xLowIndex+split;}
		}
		
	}//for (k=0; k<numberBisections; k++)
	
	
	/***************************************************************************	
	 ***** For the fit-surface case, the integrator state is stored in going 
	 ***** from the surface to the fitting point, but I am looping on j from the 
	 ***** fitting point to the surface, so:
	 ***** 
	 ***** j+1 = point closer to surface: r_old
	 ***** j   = point farther from surface: 
	 *****	   --> = r = final point storing the coefficients for the interpolation
	 ***** 
	 ***** A similar argument holds for the interval fit-center.
	 ***** 
	 ***** Visually:
	 ***** Fit -----j----(j+1) Surf
	 *****          r....r_old
	 ***** center --(j+1)---j-- Fit
	 *****          rold....r
	 ***** 
	 ***** No matter which direction I am traveling from, I need the value of Rij 
	 ***** at a given point and the espression for theta is given by:
	 ***** 
	 ***** theta = (rStep - rOld)/hUsed
	 ***** hUsed = r - rOld
	 ***** 
	 ***************************************************************************/
 	for (j=xLowIndex; j<xHighIndex; j++)
    {
		i = j;
		if (increasingR)
		{
			/* from fit to surface */
			if (rRel_and_rkf45_state_func[j+1][0] > rRelStep){break;}
		}
		else
		{
			/* from fit to center */
			if (rRel_and_rkf45_state_func[j+1][0] < rRelStep){break;}
		}
    }
	r = j;
	rOld = j+1;
	
	hUsed = rRel_and_rkf45_state_func[r][0]-rRel_and_rkf45_state_func[rOld][0];
	
	theta = (rRelStep - rRel_and_rkf45_state_func[rOld][0])/hUsed;

	i0 = 1.0 + theta*theta*(3.0-4.0*theta);
	i1 = -theta*(theta-1.0);
	i6 = -4.0*theta*theta*(theta-1.0);
	iend = theta*theta*(4.0*theta - 3.0);
	
	/* Evaluate interpolants*/
	k = 0;
	for (i = 0; i < sizeRij_calc_func; i++, k = k+2)
		Rij_from_rkf45_state[i] = i0 * Rij_calc_func[rOld][i] + 
		iend * Rij_calc_func[r][i] + 
		hUsed * i1 * rRel_and_rkf45_state_func[r][1 + k] + 
		hUsed * i6 * rRel_and_rkf45_state_func[r][2 + k];	
	
	/*fill in the matrix R*/
	k = 0;
	for (i = 0; i < rowsRic_func; i++)
		for (j = 0; j < rowsRic_func; j++, k++)
		{		
			gsl_matrix_set(matR_func,i,j, Rij_from_rkf45_state[k]);					
		}	
	
}
