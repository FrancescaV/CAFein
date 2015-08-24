#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <vector>
#include "matrix_operations.h"
#include "IOfiles.h"

using namespace std;
/*************************************************************************************
 ************ 
 ************	The routines described below handle the operations between matrices
 ************	and vectors. ALL THE MATRICES HAVE TO BEEN SQUARED!
 ************ 
 ************ 
 **************************************************************************************/



/*************************************************************************************
 *** 
 *** This routine performs the rows x columns multiplication between two matrices
 *** 
 *** Input: 
 *** - size_func = # of rows in the matrix
 *** - M1_func = matrix M1 of size (size_func, size_func)
 *** - M2_func = matrix M2 of size (size_func, size_func)
 *** - Vcoeffs_func = vector of size (size_func x size_func) to store the coefficients
 *** of the multiplication
 *** 
 *** output:
 *** - P_func = matrix product
 *** 
 *** First I calculate all the coefficients of the matrix product. 
 *** I then store them in Vcoeffs_func, and finally I fill in the matrix product.
 *** This way, I can owevwrite a matrix involved in the product without screwing up the result.
 *** 
 **************************************************************************************/
void matrix_mul_rows_by_cols(gsl_matrix *M1_func, gsl_matrix *M2_func, gsl_matrix *P_func, 
							 int size_func, gsl_vector *Vcoeffs_func)
{	
	double pij = 0.0;
	int i = 0, j = 0, k = 0, counter = 0;
	gsl_vector_set_zero(Vcoeffs_func);
	
	for(i = 0; i<size_func; i++)
		for (j = 0; j<size_func; j++)
		{
			pij = 0.0;
			for (k=0; k<size_func; k++)
			{
				pij = pij + gsl_matrix_get(M1_func,i,k)*gsl_matrix_get(M2_func,k,j);
				gsl_vector_set(Vcoeffs_func, counter, pij);
			}
			counter++;
		}

	counter = 0;
	for( i = 0; i<size_func; i++)
		for ( j = 0; j<size_func; j++)
		{
			gsl_matrix_set(P_func,i,j, gsl_vector_get(Vcoeffs_func, counter));
			counter++;
		}

}//matrix_mul_rows_by_cols

/*************************************************************************************
 *** 
 *** This routine performs the sum of two matrices
 *** 
 *** Input: 
 *** - size_func = # of rows in the matrix
 *** - M1_func = matrix M1 of size (size_func, size_func)
 *** - M2_func = matrix M2 of size (size_func, size_func)
 *** 
 *** output:
 *** - S_func = matrix sum
 *** 
 **************************************************************************************/
void matrix_sum (gsl_matrix *M1_func, gsl_matrix *M2_func, gsl_matrix *S_func, int size_func)
{
	int i=0, j=0;

	for (i = 0; i<size_func; i++)
		for(j = 0; j<size_func; j++)
			gsl_matrix_set(S_func,i,j, gsl_matrix_get(M1_func,i,j) + gsl_matrix_get(M2_func,i,j));

}//matrix_sum
/*************************************************************************************
 *** 
 *** This routine performs the difference of two matrices
 *** 
 *** Input: 
 *** - size_func = # of rows in the matrix
 *** - M1_func = matrix M1 of size (size_func, size_func)
 *** - M2_func = matrix M2 of size (size_func, size_func)
 *** 
 *** output:
 *** - D_func = matrix difference
 *** 
 **************************************************************************************/
void matrix_diff (gsl_matrix *M1_func, gsl_matrix *M2_func, gsl_matrix *D_func, int size_func)
{
	int i=0, j=0;
	
	for (i = 0; i<size_func; i++)
		for(j = 0; j<size_func; j++)
			gsl_matrix_set(D_func,i,j, gsl_matrix_get(M1_func,i,j) - gsl_matrix_get(M2_func,i,j));
}//matrix_diff
/*************************************************************************************
 *** 
 *** This routine extracts from a matrix of size (sizeMi_func*2, sizeMi_func*2)
 *** 4 submatrices of size (sizeMi_func, sizeMi_func)
 *** 
 *** Input: 
 *** - sizeMi_func = # of rows for submatrices
 *** - M_func = total matrix
 *** 
 *** output:
 *** - M1_func = submatrix 1
 *** - M2_func = submatrix 2
 *** - M3_func = submatrix 3
 *** - M4_func = submatrix 4
 *** 
 *** such that:
 *** 
 ***	|M1  M2|
 *** M =|	   |
 ***	|M3  M4|
 *** 
 *** 
 **************************************************************************************/
void extract_submatrices(gsl_matrix * M_func, gsl_matrix * M1_func, gsl_matrix * M2_func, 
						 gsl_matrix * M3_func,gsl_matrix * M4_func, int sizeMi_func)
{	
	for (int k=0; k<sizeMi_func; k++)
		for (int j=0; j<sizeMi_func; j++)
		{
			gsl_matrix_set(M1_func,k,j, gsl_matrix_get(M_func,k,j));
			gsl_matrix_set(M2_func,k,j, gsl_matrix_get(M_func,k,j+sizeMi_func));
			gsl_matrix_set(M3_func,k,j, gsl_matrix_get(M_func,k+sizeMi_func,j));
			gsl_matrix_set(M4_func,k,j, gsl_matrix_get(M_func,k+sizeMi_func,j+sizeMi_func));
			
		}	

}//extract_submatrices
/*************************************************************************************
 *** 
 *** This routine performs the rows x columns multiplication between a matrix and a vector
 *** 
 *** Input: 
 *** - size_func = # of rows in the matrix = size of vector
 *** - M_func = matrix of size (size_func, size_func)
 *** - V_func = vector of size size_func
 *** - Vcoeffs_func = vector of size size_func to store the coefficients
 *** output:
 *** - P_func = vector product of size size_func
 *** 
 *** First I calculate all the coefficients of the product. 
 *** I then store them in Vcoeffs_func, and finally I fill in the vector product.
 *** This way, I can owevwrite a vector involved in the product without screwing up the result. 
 **************************************************************************************/
void mul_matrix_by_vec(gsl_matrix *M_func, gsl_vector *V_func, gsl_vector *P_func, gsl_vector *Vcoeffs_func,
					   int size_func)
{
	int i = 0, j = 0;
	double p_i = 0.0;
	
	for(i = 0; i<size_func; i++)
	{
		p_i = 0.0;
		for(j = 0; j<size_func; j++)
		{
			p_i = p_i + gsl_matrix_get(M_func,i,j)*gsl_vector_get(V_func,j);
			gsl_vector_set(Vcoeffs_func, i, p_i);
		}
	}
	
	for (i = 0; i<size_func; i++)
		gsl_vector_set(P_func,i, gsl_vector_get(Vcoeffs_func,i));
	
		
}//mul_matrix_by_vec
/*************************************************************************************
 *** 
 *** This routine calculates the determinant of a matrix
 *** 
 *** Input: 
 *** - size_func = # of rows in the matrix
 *** - M_func = matrix of size (size_func, size_func)
 *** - p = gsl_permutation used to calculate the determinant
 *** - DummyMatrix = storage matrix of size (size_func, size_func)
 *** 
 *** output:
 *** - determinant = determinant of M
 *** 
 **************************************************************************************/
double determinantNbyN(gsl_matrix *M_func, int size_func, gsl_permutation *p, gsl_matrix *DummyMatrix)
{
	double determinant = 0.0;
	
	
	if (size_func == 2)
	{
		determinant = gsl_matrix_get(M_func,0,0)*gsl_matrix_get(M_func,1,1) - 
		gsl_matrix_get(M_func,0,1)*gsl_matrix_get(M_func,1,0);
	}
	else if (size_func >2)
	{
		int s;
		gsl_permutation_init(p);

		gsl_matrix_memcpy(DummyMatrix, M_func);
		gsl_linalg_LU_decomp (M_func, p, &s);
		determinant = gsl_linalg_LU_det(M_func, s);
		gsl_matrix_memcpy(M_func, DummyMatrix);	
	}
	return determinant;
	
}//determinantNbyN
/*************************************************************************************
 *** 
 *** This routine calculates the determinant of the Riccati matrix R in its complex version
 *** 
 *** R = Rrr+i Rir
 *** 
 *** and R is 3x3 (non adiabatic case). 
 *** The total determinant has a real and imaginary part.
 *** 
 *** Input: 
 *** - Rrr_func = real matrix R of size (3,3)
 *** - Rir_func = imaginary matrix R of size (3,3)
 *** 
 *** output:
 *** - vecDetDiff_r_i_func = vector of size 2 that stores the real and imaginary part of the
 *** determinant (elements #0 and #1, respectively)
 **************************************************************************************/
void complexDeterminant3by3(gsl_matrix *Rrr_func, gsl_matrix *Rir_func, gsl_vector *vecDetDiff_r_i_func)
{
	
	double Rr00 = 0.0, Rr01 = 0.0, Rr02 = 0.0, Rr10 = 0.0, Rr11 = 0.0, Rr12 = 0.0, Rr20 = 0.0, Rr21 = 0.0, Rr22 = 0.0, 
	Ri00 = 0.0, Ri01 = 0.0, Ri02 = 0.0, Ri10 = 0.0, Ri11 = 0.0, Ri12 = 0.0, Ri20 = 0.0, Ri21 = 0.0, Ri22 = 0.0;
	
	Rr00 = gsl_matrix_get(Rrr_func,0,0);
	Rr01 = gsl_matrix_get(Rrr_func,0,1);
	Rr02 = gsl_matrix_get(Rrr_func,0,2);
	Rr10 = gsl_matrix_get(Rrr_func,1,0);
	Rr11 = gsl_matrix_get(Rrr_func,1,1);
	Rr12 = gsl_matrix_get(Rrr_func,1,2);
	Rr20 = gsl_matrix_get(Rrr_func,2,0);
	Rr21 = gsl_matrix_get(Rrr_func,2,1);
	Rr22 = gsl_matrix_get(Rrr_func,2,2);
	
	Ri00 = gsl_matrix_get(Rir_func,0,0);
	Ri01 = gsl_matrix_get(Rir_func,0,1);
	Ri02 = gsl_matrix_get(Rir_func,0,2);
	Ri10 = gsl_matrix_get(Rir_func,1,0);
	Ri11 = gsl_matrix_get(Rir_func,1,1);
	Ri12 = gsl_matrix_get(Rir_func,1,2);
	Ri20 = gsl_matrix_get(Rir_func,2,0);
	Ri21 = gsl_matrix_get(Rir_func,2,1);
	Ri22 = gsl_matrix_get(Rir_func,2,2);
	
	/*real part of the determinant*/
	gsl_vector_set(vecDetDiff_r_i_func,0, -(Ri22*(Ri11*Rr00 - Ri10*Rr01 - Ri01*Rr10 + Ri00*Rr11)) +
				   Ri21*(Ri12*Rr00 - Ri10*Rr02 - Ri02*Rr10 + Ri00*Rr12) - 
				   Ri20*(Ri12*Rr01 - Ri11*Rr02 - Ri02*Rr11 + Ri01*Rr12) + 
				   (Ri02*Ri11 - Ri01*Ri12 - Rr02*Rr11 + Rr01*Rr12)*Rr20 - 
				   (Ri02*Ri10 - Ri00*Ri12 - Rr02*Rr10 + Rr00*Rr12)*Rr21 + 
				   (Ri01*Ri10 - Ri00*Ri11 - Rr01*Rr10 + Rr00*Rr11)*Rr22);
		
	
	/*imaginary part of the determinant*/
	gsl_vector_set(vecDetDiff_r_i_func,1, Ri22*(Ri01*Ri10 - Ri00*Ri11 - Rr01*Rr10 + Rr00*Rr11) - 
				   Ri21*(Ri02*Ri10 - Ri00*Ri12 - Rr02*Rr10 + Rr00*Rr12) + 
				   Ri20*(Ri02*Ri11 - Ri01*Ri12 - Rr02*Rr11 + Rr01*Rr12) + 
				   (Ri12*Rr01 - Ri11*Rr02 - Ri02*Rr11 + Ri01*Rr12)*Rr20 - 
				   (Ri12*Rr00 - Ri10*Rr02 - Ri02*Rr10 + Ri00*Rr12)*Rr21 + 
				   (Ri11*Rr00 - Ri10*Rr01 - Ri01*Rr10 + Ri00*Rr11)*Rr22);
	
}//complexDeterminant3by3
/*************************************************************************************
 *** 
 *** This routine calculates the inverse of a matrix
 *** 
 *** Input: 
 *** - size_func = # of rows in the matrix
 *** - M_func = matrix of size (size_func, size_func)
 *** - p = gsl_permutation used to calculate the inverse
 *** - DummyMatrix = storage matrix of size (size_func, size_func)
 *** 
 *** output:
 *** - M_inv_func = inverse of M
 *** 
 **************************************************************************************/

void InverseMatrix(gsl_matrix *M_func, gsl_matrix *M_inv_func, int size_func, gsl_permutation *p, gsl_matrix *DummyMatrix)
{
	
	int s = 0;
	double detM = 0.0;
	
	if (size_func == 2)
	{
	
		double m00 = 0.0, m01 = 0.0, m10 = 0.0, m11 = 0.0, detM = 0.0;
		
		m00 = gsl_matrix_get(M_func,0,0);
		m01 = gsl_matrix_get(M_func,0,1);
		m10 = gsl_matrix_get(M_func,1,0);
		m11 = gsl_matrix_get(M_func,1,1);
		
		detM = m00*m11 - m01*m10;
		
		if (detM == 0.0)
			errorMessageAndExit("matrix_operations.c (InverseMatrix, size 2)", " SINGULAR MATRIX ");
		
		gsl_matrix_set(M_inv_func,0,0, m11/detM);
		gsl_matrix_set(M_inv_func,0,1, -m01/detM);
		gsl_matrix_set(M_inv_func,1,0, -m10/detM);
		gsl_matrix_set(M_inv_func,1,1, m00/detM);
	}
	else if(size_func == 3)
	{
	
		double m00 = 0.0, m01 = 0.0, m02 = 0.0, 
		m10 = 0.0, m11 = 0.0, m12 = 0.0,
		m20 = 0.0, m21 = 0.0, m22 = 0.0;
		
		m00 = gsl_matrix_get(M_func,0,0);
		m01 = gsl_matrix_get(M_func,0,1);
		m02 = gsl_matrix_get(M_func,0,2);
		
		m10 = gsl_matrix_get(M_func,1,0);
		m11 = gsl_matrix_get(M_func,1,1);
		m12 = gsl_matrix_get(M_func,1,2);
		
		m20 = gsl_matrix_get(M_func,2,0);
		m21 = gsl_matrix_get(M_func,2,1);
		m22 = gsl_matrix_get(M_func,2,2);
		
		detM = -(m02*m11*m20) + m01*m12*m20 + m02*m10*m21 - m00*m12*m21 - m01*m10*m22 + m00*m11*m22;
		
		if (detM == 0.0)
			errorMessageAndExit("matrix_operations.c (InverseMatrix, size 3)", " SINGULAR MATRIX ");

		gsl_matrix_set(M_inv_func,0,0, (-(m12*m21) + m11*m22)/detM);
		gsl_matrix_set(M_inv_func,0,1, (m02*m21 - m01*m22)/detM);
		gsl_matrix_set(M_inv_func,0,2, (-(m02*m11) + m01*m12)/detM);
		
		gsl_matrix_set(M_inv_func,1,0, (m12*m20 - m10*m22)/detM);
		gsl_matrix_set(M_inv_func,1,1, (-(m02*m20) + m00*m22)/detM);
		gsl_matrix_set(M_inv_func,1,2, (m02*m10 - m00*m12)/detM);
		
		gsl_matrix_set(M_inv_func,2,0, (-(m11*m20) + m10*m21)/detM);
		gsl_matrix_set(M_inv_func,2,1, (m01*m20 - m00*m21)/detM);
		gsl_matrix_set(M_inv_func,2,2, (-(m01*m10) + m00*m11)/detM);	
	}
	else if(size_func >3)
	{		
		gsl_matrix_memcpy(DummyMatrix, M_func);
		gsl_linalg_LU_decomp (M_func, p, &s);    
		gsl_linalg_LU_invert (M_func, p, M_inv_func);
		gsl_matrix_memcpy(M_func, DummyMatrix);
	}	
	
}//InverseMatrix

/*************************************************************************************
 *** 
 *** This routine calculates the sum of two vectors:
 *** 
 *** Input: 
 *** - size_func = size of the vector
 *** - V1_func = vector 1
 *** - V2_func = vector 2
 *** 
 *** output:
 *** - S_func = vector sum V1 + V2.
 *** 
 **************************************************************************************/
void vector_sum (gsl_vector *V1_func, gsl_vector *V2_func, gsl_vector *S_func, int size_func)
{
	int i=0;
	
	for (i = 0; i<size_func; i++)
		gsl_vector_set(S_func, i, gsl_vector_get(V1_func,i) + gsl_vector_get(V2_func,i));
	
}//vector_sum
