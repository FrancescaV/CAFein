#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <limits.h>
#include <gsl/gsl_matrix.h>
#include "dma.h"
#include "SteffenInterp.h"

using namespace std;

double signum(double n)
{
	double sign = 0.0;
	
	if (n >= 0.0) {sign = 1.0;}
	else {sign = -1.0;}

	return sign;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This is stupid. I couldn't compile using min and max algorithm...so i had to create one!
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double findMin(double a, double b)
{
	double min = 0.0;
	if (a<b) {min = a;}
	if (a>b) {min = b;}
	return min;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void calcSteffenInterp(int N,double *x, double *y,double **c)
{
/*
!**** Calculate a piecewise cubic interpolant of a tabulated function (x_i,y_i) that is monotonic on
!**** each interval (x_i, x_i+1) for i=0,2,...,N-2. On each interval, the interpolant is of the form
!**** 
!****    f_i(x) = c_(1,i)*(x-x_i)^3 + c_(2,i)*(x-x_i)^2 
!****           + c_(3,i)*(x-x_i) + c_(4,i)
!**** 
!**** The interpolant and its first 3 derivatives can be evaluated by calling the SUBROUTINE 
!**** EvalSteffenInterp. Note however that the interpolation scheme only guarantees continuity
!**** of the function and its first derivative, not of the higher order derivatives. 
!**** The method is described in detail in 
!****    Steffen, M., 1990, A&A 239, 443
!**** 
!**** Input parameters
!**** ----------------
!****    N: The number of data points
!****    x: 1-D array of size N containing the independent variable
!****    y: 1-D array of size N containing the dependent variable
!**** 
!**** Output parameters
!**** -----------------
!****    c: 2-D array of size 4x(N-1) containing the coefficients 
!****       of the interpolant on each interval (x_i,x_i+1) for 
!****       i=0,2,...,N-2
!**** 
!**** Coded by Bart Willems on 08/28/06
!**** Rewritten in c++ by Francesca Valsecchi 10/30/08
!****
*/

	int i = 0;
    double 	p;
    double	*dy, *h, *s;
    	dy = dfun(N);
    	h = dfun(N-1);
    	s = dfun(N-1);
// Determine the sizes h_i of the interval (x_i,x_i+1) and the slope of the secants through (x_i,y_i) 
// and (x_i+1,y_i+1) for i=0,2,...,N-2 [Eqs. (6) and (7) in Steffen (1990)]
	for (i=0; i<N-1; i++)
     {
         h[i] = x[i+1]-x[i];
         s[i] = (y[i+1]-y[i])/h[i];
     }
// Determine the derivative at (x_i,y_i) for i=0,2,...,N-1 [Eqs. (11), (26), (27) in Steffen (1990)]
	
    p = s[0]*(1.0+(h[0]/(h[0]+h[1]))) - s[1]*h[0]/(h[0]+h[1]);
    dy[0] = (signum(p)+signum(s[0]))*
    		findMin(fabs(s[0]),0.5*fabs(p));
	
	for (i=1; i<N-1; i++)
	{
		p = (s[i-1]*h[i]+s[i]*h[i-1])/(h[i-1]+h[i]);
        dy[i] = (signum(s[i-1])+signum(s[i]))*
        		findMin(fabs(s[i-1]),findMin(fabs(s[i]),0.5*fabs(p)));

	}
    p = s[N-2]*(1.0+(h[N-2]/(h[N-2]+h[N-3]))) - s[N-3]*h[N-2]/(h[N-2]+h[N-3]);
	dy[N-1] = (signum(p) + signum(s[N-2]))*
			findMin(fabs(s[N-2]),0.5*fabs(p));

// Determine the coefficients of the cubic interpolants on 
// the intervals (x_i,x_i+1) for i=0,2,...,N-2
// [Eqs. (2)-(5) in Steffen (1990)]
	for(i=0; i<N-1; i++)
    {
		c[0][i] = (dy[i]+dy[i+1]-2.0*s[i])/(h[i]*h[i]);
		c[1][i] = (3.0*s[i]-2.0*dy[i]-dy[i+1])/h[i];
		c[2][i] = dy[i];			
        c[3][i] = y[i];
    }
	
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void evalSteffenInterp(int N,double *x, double **c, double x0,double *f)
{
/*
!**** Evaluate the Steffen interpolant and its first 3 derivatives at x0. On each interval 
!**** [x_i,x_i+1] with i=0,2,...,N-2, the interpolant is of the form
!**** 
!****    f_i(x) = c_(1,i)*(x-x_i)^3 + c_(2,i)*(x-x_i)^2 
!****           + c_(3,i)*(x-x_i) + c_(4,i)
!**** 
!**** The coefficients c_(1,i), c_(2,i), c_(3,i), c_(4,i) are determined by calling the SUBROUTINE 
!**** calcSteffenInterpol. Note however that the interpolation scheme only guarantees continuity 
!**** of the function and its first derivative, not of the higher order derivatives. The method is 
!**** described in detail in
!**** 
!****    Steffen, M., 1990, A&A 239, 443
!**** 
!**** Input parameters
!**** ----------------
!****    N: The number of data points
!****    x: 1-D array of size N containing the independent variable
!****    c: 2-D array of size 4x(N-1) containing the coefficients 
!****       of the interpolant on each interval (x_i,x_i+1) with
!****       i=0,2,...,N-2
!****    x0: The point at which the interpolant and its derivatives
!****        is to be evaluated
!**** 
!**** Output parameters
!**** -----------------
!****    f: 1-D array of size 4 containing the value of the interpolant
!****       and its first 3 derivatives at x0 
!****
!****       -> f(0) contains the function value
!****       -> f(1) contains the 1st derivative
!****       -> f(2) contains the 2nd derivative
!****       -> f(3) contains the 3rd derivative
!**** 
!**** Coded by Bart Willems on 08/28/06
!**** Rewritten in c++ by Francesca Valsecchi 10/30/08
!****
*/
	int i = 0, k = 0;
// Check if x0 is in the range of the (x_i) data set
	if (x0 < x[0] || x0 > x[N-1])
	{
		cout << "FATAL ERROR: Extrapolation not allowed in " 
			 << "SUBROUTINE EvalSteffenInterp!'" << endl;	         	
	   	exit(EXIT_FAILURE);
    }
// Locate x0 in the (x_i) data set
	int xHighIndex = N-1, xLowIndex = 0;
	int difference = 0, split = 0;
	int numberBisections = 100;
	
	for (k=0; k<numberBisections; k++)
	{

		/* determining the interval where the crossing happens by bisection*/
		difference = xHighIndex-xLowIndex;
		if (difference < 10) {break;}					

		if (difference%2 == 0) {;}
		else {difference = difference +1;}
		
		split = difference/2;

		if (x0 >= x[xLowIndex] && x0 < x[xLowIndex+split]) {xHighIndex = xLowIndex+split; }
		else {xLowIndex = xLowIndex+split;}

	}//for (k=0; k<numberBisections; k++)

	for (int j=xLowIndex; j<xHighIndex; j++)
    {
         i = j;
         if (x0 < x[j+1]) {break;}         
    }

// Evaluate the interpolant and its first 3 derivatives at x0
      f[0] = c[0][i]*pow((x0-x[i]),3.0) + c[1][i]*pow((x0-x[i]), 2.0) + 
      		 c[2][i]*(x0-x[i]) + c[3][i];
      f[1] = 3.0*c[0][i]*pow((x0-x[i]),2.0) + 2.0*c[1][i]*(x0-x[i]) + c[2][i];
      f[2] = 6.0*c[0][i]*(x0-x[i]) + 2.0*c[1][i];
      f[3] = 6.0*c[0][i];


}

void calcSteffenInterp_RiccatiElmnts(int Ny, int position_xIn, int position_xFin, 
									 double **dy, double **h, double **s, double *x, double **y, double **c)
{
	/*
	 !**** Calculate a piecewise cubic interpolant of tabulated functions across a specific interval (xIn,xFin)
	 !**** that is monotonic on each interval (x_i, x_i+1) for i=positionXin........positionXfin-2.
	 !**** On each interval, the interpolant is of the form
	 !****
	 !****    f_i(x) = c_(1,i)*(x-x_i)^3 + c_(2,i)*(x-x_i)^2 
	 !****           + c_(3,i)*(x-x_i) + c_(4,i)
	 !**** 
	 !**** The interpolant and its first derivative can be evaluated by calling the SUBROUTINE 
	 !**** evalSteffenInterp_RiccatiElmnts. The interpolation scheme guarantees continuity
	 !**** of the function and its first derivative. 
	 !**** The method is described in detail in 
	 !****    Steffen, M., 1990, A&A 239, 443
	 !**** 
	 !**** Input parameters
	 !**** ----------------
	 !****    Ny: The number of functions to fit
	 !****    position_xIn: position of beginning interpolation interval
	 !****    position_xFin: position of end interpolation interval
	 !****    dy: 2-D array of size (Nx, Ny)used by Steffen to perform calculation
	 !****    h: 2-D array of size (Nx-1, Ny)used by Steffen to perform calculation
	 !****    s: 2-D array of size (Nx-1, Ny)used by Steffen to perform calculation
	 !****    x: 1-D array of size Nx containing the independent variable
	 !****    y: 2-D array of size (Nx, Ny) containing the dependent variable
	 !**** 
	 !**** Output parameters
	 !**** -----------------
	 !****    c: 2-D array of size ( 4*Ny, Nx-1) containing the coefficients 
	 !****       of the interpolant on each interval (x_i,x_i+1) for 
	 !****       i=positionXin........positionXfin-2
	 !**** 
	 !**** Written in c++ by Francesca Valsecchi 03/16/12
	 !****
	 */
	
    double 	p;
	int i = 0, j = 0, k = 0, 
	N = position_xFin - position_xIn + 1;
	/*****************************************************************************************
	 *******  Determine the sizes h_i of the interval (x_i,x_i+1) and the slope of the
	 *******  secants through (x_i,y_i) and (x_i+1,y_i+1) for i=0...,N-2
	 *******  [Eqs. (6) and (7) in Steffen (1990)]
	 ******************************************************************************************/	 
	for (i = 0; i < N-1; i++)
		for (j = 0; j < Ny; j++)
		{
			k = position_xIn + i;
			
			h[k][j] = x[k+1]-x[k];
			s[k][j] = (y[k+1][j]-y[k][j])/h[k][j];
		}
	
	/***********************************************************************
	 *******  Determine the derivative at (x_i,y_i) for i=1,2,...,N-1
	 *******  [Eqs. (11), (26), (27) in Steffen (1990)]
	************************************************************************/
	/*first row*/
	for (j = 0; j < Ny; j++)
	{
		k = position_xIn;

		p = s[k][j]*(1.0+(h[k][j]/(h[k][j]+h[k+1][j]))) - s[k+1][j]*h[k][j]/(h[k][j]+h[k+1][j]);
		dy[k][j] = (signum(p)+signum(s[k][j]))*findMin(fabs(s[k][j]),0.5*fabs(p));
	}
	/*remaining elements*/
	for (i = 1; i < N-1; i++)
		for (j = 0; j < Ny; j++)
		{			
			k = position_xIn + i;

			p = (s[k-1][j]*h[k][j]+s[k][j]*h[k-1][j])/(h[k-1][j]+h[k][j]);
			dy[k][j] = (signum(s[k-1][j])+signum(s[k][j]))*
			findMin(fabs(s[k-1][j]),findMin(fabs(s[k][j]),0.5*fabs(p)));
		}

	/*last element*/
	for (j = 0; j < Ny; j++)
	{
		k = position_xIn + N;
		
		p = s[k-2][j]*(1.0+(h[k-2][j]/(h[k-2][j]+h[k-3][j]))) - s[k-3][j]*h[k-2][j]/(h[k-2][j]+h[k-3][j]);
		dy[k-1][j] = (signum(p) + signum(s[k-2][j]))*findMin(fabs(s[k-2][j]),0.5*fabs(p));
	}	
	/***********************************************************************
	 *******  Determine the coefficients of the cubic interpolants on the intervals
	 *******  (x_i,x_i+1) for i=0...,N-2
	 *******  [Eqs. (2)-(5) in Steffen (1990)]
	 ************************************************************************/
	for(i = 0; i <  N-1; i++)
		for (j = 0; j < Ny; j++)
		{
			k = position_xIn + i;
			
			c[0 + 4*j][k] = (dy[k][j]+dy[k+1][j]-2.0*s[k][j])/(h[k][j]*h[k][j]);
			c[1 + 4*j][k] = (3.0*s[k][j]-2.0*dy[k][j]-dy[k+1][j])/h[k][j];
			c[2 + 4*j][k] = dy[k][j];				
			c[3 + 4*j][k] = y[k][j];
		}
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void evalSteffenInterp_RiccatiElmnts(int Ny, int position_xIn, int position_xFin, 
									 double *x, double **c, double x0,double *f, gsl_matrix* matY, int rowsY)
{
	/*
	 !**** Evaluate the Steffen interpolant and its first derivative at x0. On each interval 
	 !**** [x_i,x_i+1] with i=position_xIn,...,position_xFin-2, the interpolant is of the form
	 !**** 
	 !****    f_i(x) = c_(1,i)*(x-x_i)^3 + c_(2,i)*(x-x_i)^2 
	 !****           + c_(3,i)*(x-x_i) + c_(4,i)
	 !**** 
	 !**** The coefficients c_(1,i), c_(2,i), c_(3,i), c_(4,i) are determined by calling the SUBROUTINE 
	 !**** calcSteffenInterp_RiccatiElmnts. The interpolation scheme guarantees continuity 
	 !**** of the function and its first derivative. The method is 
	 !**** described in detail in
	 !**** 
	 !****    Steffen, M., 1990, A&A 239, 443
	 !**** 
	 !**** Input parameters
	 !**** ----------------
	 !****    Ny: The number of functions to fit
	 !****    position_xIn: position of beginning interpolation interval
	 !****    position_xFin: position of end interpolation interval
	 !****    x: 1-D array of size Nx containing the independent variable
	 !****    c: 2-D array of size ( 4*Ny, Nx-1) containing the coefficients 
	 !****       of the interpolant on each interval (x_i,x_i+1) for 
	 !****       i=positionXin........positionXfin-2
	 !****    x0: The point at which the interpolant and its derivatives
	 !****        is to be evaluated
	 !****    rowsY: rows of matrix matY
	 !**** 
	 !**** Output parameters
	 !**** -----------------
	 !****    f: 1-D array of size 2*Ny containing the value of the interpolant
	 !****       and its first derivative at x0 
	 !****    matY: the riccati matrix containing Rij (or y variables) 
	 !****
	 !****       -> f(0) contains the function value
	 !****       -> f(1) contains the 1st derivative
	 !**** 
	 !**** Written in c++ by Francesca Valsecchi 03/16/12
	 !****
	 */
	int i = 0, k = 0, j = 0, increasingX = 0, halfY = rowsY/2;;
	
	if (x[position_xFin] > x[position_xIn]){increasingX = 1;}
	
	/* Check if x0 is in the range of the (x_i) data set*/
	
	if (increasingX)
	{
		if (x0 < x[position_xIn] || x0 > x[position_xFin])
		{
			cout << "FATAL ERROR: Extrapolation not allowed in " << "SUBROUTINE EvalSteffenInterp!'" << endl;	         	
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		if (x0 > x[position_xIn] || x0 < x[position_xFin])
		{
			cout << "FATAL ERROR: Extrapolation not allowed in " << "SUBROUTINE EvalSteffenInterp!'" << endl;	         	
			exit(EXIT_FAILURE);
		}
	}
	
	/* Locate x0 in the (x_i) data set*/
	int xHighIndex = position_xFin, xLowIndex = position_xIn;
	int difference = 0, split = 0;
	
	int numberBisections = 10000;
	
	for (k=0; k<numberBisections; k++)
	{
		/* determining the interval where the crossing happens by bisection*/
		difference = xHighIndex-xLowIndex;
		if (difference < 10) {break;}					
		
		if (difference%2 == 0) {;}
		else {difference = difference +1;}
		
		split = difference/2;
		
		if(increasingX)
		{
			if (x0 >= x[xLowIndex] && x0 < x[xLowIndex+split]) {xHighIndex = xLowIndex+split; }
			else {xLowIndex = xLowIndex+split;}
		}			
		else
		{
			if (x0 <= x[xLowIndex] && x0 > x[xLowIndex+split]) {xHighIndex = xLowIndex+split; }
			else {xLowIndex = xLowIndex+split;}
		}

	}//for (k=0; k<numberBisections; k++)
	
	for (j=xLowIndex; j<xHighIndex; j++)
    {
		i = j;
		if (increasingX)
		{
			if (x[j+1] > x0){break;}
		}
		else
		{
			if (x[j+1] < x0){break;}
		}
    }
	
	// Evaluate the interpolant and its first derivative at x0
	for (j = 0; j < Ny; j++)
	{
		f[0+2*j] = c[0+4*j][i]*pow((x0-x[i]),3.0) + c[1+4*j][i]*pow((x0-x[i]), 2.0) + c[2+4*j][i]*(x0-x[i]) + c[3+4*j][i];
		f[1+2*j] = 3.0*c[0+4*j][i]*pow((x0-x[i]),2.0) + 2.0*c[1+4*j][i]*(x0-x[i]) + c[2+4*j][i];
		
	}	
	/*filling in R both in adiabatic and non adiabatic case*/
	k = 0;
	if(Ny == 4)
	{
		for (i = 0; i < Ny/2; i++)
			for (j = 0; j < Ny/2; j++, k = k+2)
				gsl_matrix_set(matY,i,j, f[k]);		
	}
	else
	{
		for (i = 0; i < halfY; i++)
			for (j = 0; j < rowsY; j++, k = k+2)
			{		
				gsl_matrix_set(matY,i,j, f[k]);		
				
				if(j < halfY){gsl_matrix_set(matY,i+halfY,j+halfY, f[k]);}
				else{gsl_matrix_set(matY,i+halfY,j-halfY, -f[k]);}
			}
	}
	
}
