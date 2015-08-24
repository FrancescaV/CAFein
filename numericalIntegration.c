#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include "numericalIntegration.h"
#include "IOfiles.h"
using namespace std;

/*************************************************************************************
 ************ 
 ************	The routines described below are used to perform numerical integration
 ************ 
 **************************************************************************************/


/*************************************************************************************
 ****** This routine is based on the algorithm presented by Gill & miller
 ****** 
 ****** http://comjnl.oxfordjournals.org/content/15/1/80.full.pdf
 ****** 
 ****** This subroutine evaluates the integral from x(ia) to x(ib) of a function
 ****** whose values at points x(i), i=0(1)n, are stored in f(i), i=0(1)n. The
 ****** points x(i), which need not be equally spaced, are assumed to be distinct
 ****** and arranged in ascending or descending order. The integration is performed
 ****** using a 4-point rule over each interval. If any of the relations n < 3,
 ****** ib > n, ia < 0, ia > ib are true, the subroutine fails. If n > 3, an
 ****** estimation of the error is assigned to er. If n = 3, no indication of the
 ****** error is obtainable, and er=0.
 ****** 
 ****** GILL, P.E. and MILLER, G.F.
 ****** An algorithm for the integration of unequally spaced data.
 ****** Computer Journal 15, pp. 80-83, 1972.
 ****** 
 ******  
 ****** Note:
 ****** ------
 ****** In their paper, Gill and Miller do not add the quantity 'er' to 'intt'
 ****** before return. However, extensive tests have shown that a dramatic reduction
 ****** in the error often results from such an addition. In other cases, it does
 ****** not make an improvement, but these tend to be cases of low accuracy in
 ****** which the modified answer is not significantly inferior to the unmodified
 ****** one. The user has the option of recovering the Gill-Miller answer by
 ****** subtracting 'er' from 'intt' on return from the routine.
 ****** 
 ****** 
 ******	- N: number of data points
 ******	- ia: Starting index
 ******	- ib: End index
 ****** - err: error
 ****** 
 **************************************************************************************/

double FourPt_DoubleStar(double *x, double *f, int N, int ia, int ib, double &error)
{
	
	double integral = 0.0;
	/*Termination conditions*/
	if(N < 3 || ia > ib || ia < 0 ||ib > N)
	{
		integral = 0.0;
		error = 0.0;		
		errorMessageAndExit("numericalIntegration.c", "Bad input in FourPt_DoubleStar");
	}
	else if (ia == ib) 
	{
		cout << "ia == ib ==> error and integral = 0.0" << endl;
		integral = 0.0;
		error = 0.0;
	}
	else
	{
		double e = 0.0, h1 = 0.0, h2 = 0.0, h3 = 0.0, h4 = 0.0, 
		r1 = 0.0, r2 = 0.0, r3 = 0.0, r4 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0, 
		c = 0.0, s = 0.0, intP = 0.0;

		int i = 0, j = 0, k = 0;
		
		if (ia == (N-1) && N > 3){j = N - 2;}
		else if (ia > 2){j = ia;}
		else {j = 2;}
	
		if (ib == 1 && N > 3){k = 3;}
		else if (N > ib + 2){k = ib + 1;}
		else {k = N - 1;}
		
		for (i = j; i <= k; i++)
		{
			if (i == j)
			{
			
				h2 = x[j-1] - x[j-2];
				d3 = (f[j-1] - f[j-2])/h2;
				
				h3 = x[j] - x[j-1];
				d1 = (f[j] - f[j-1])/h3;
				
				h1 = h2 + h3;
				d2 = (d1 - d3)/h1;
				
				h4 = x[j + 1] - x[j];
				r1 = (f[j + 1] - f[j])/h4;
				
				r2 = (r1 - d1)/(h4 + h3);
				
				h1 = h1 + h4;
				
				r3 = (r2 - d2)/h1;
				
				if (ia == 0)
				{
					intP = h2*(f[0]+h2*(d3/2.0-h2*(d2/6.0-(h2 + 2.0*h3)*r3/12.0)));
					s = -pow(h2, 3.0)*(h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0;
				}
			}//if (i == j)
			else
			{
				h4 = x[i + 1] - x[i];
				r1 = (f[i + 1] - f[i])/h4;
				r4 = h4 + h3;
				r2 = (r1 - d1)/r4;
				r4 = r4 + h2;
				r3 = (r2 - d2)/r4;
				r4 = r4 + h1;
				r4 = (r3 - d3)/r4;
			}
			
			if( i <= ib && i > ia)
			{
				intP = intP+h3*((f[i]+f[i-1])/2.0-h3*h3*(d2+r2+(h2-h4)*r3)/12.0);
				c = pow(h3,3.0)*(2.0*h3*h3+5.0*(h3*(h4+h2)+2.0*h4*h2))/120.0;
				e = e+(c+s)*r4;
				
				if(i == j){s = 2.0*c+s;}
				else{s = c;}
			}//if( i <= ib && i > ia)
			else{e = e+r4*s;}
			
			if (i == k)
			{
				if(ib == N)
				{
					intP = intP+h4*(f[N]-h4*(r1/2.0+h4*(r2/6.0+(2.0*h3+h4)*r3/12.0)));
					e = e-pow(h4, 3.0)*r4*(h4*(3.0*h4+5.0*h2)+10.0*h3*(h2+h3+h4))/60.0;					
				}//if(ib == N)
				if(ib >= N - 1){e = e+s*r4;}				
			}//if (i == k)
			else
			{
				h1 = h2;
				h2 = h3;
				h3 = h4;
				d1 = r1;
				d2 = r2;
				d3 = r3;			
			}
		}//for (i = j; i <= k; i++)
		error = e;
		integral = intP+error;
	}
	return integral;

}

