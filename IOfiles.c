#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ostream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "IOfiles.h"
using namespace std;

/*************************************************************************************
 ************ 
 ************	The routines described below handle the input and output files.
 ************ 
 **************************************************************************************/

/*************************************************************************************
 ************ 
 ************	This routine takes an input file and stores the elements of the file
 ************	in a matrix, outputting also the number of lines in the file.
 ************	
 ************	Input:
 ************	- fileName = address of the input file
 ************	Output:
 ************	- table_func = matrix containing the input file elements
 ************	- num_lines = number of lines in the file
 ************	
 **************************************************************************************/
int readInputFile(const char *fileName, vector<vector<double> > &table_func)
{
	ifstream file;
	file.open(fileName); 	
	
	int num_lines = 0;
	double data = 0.0;
	
	while(file)
	{
		string eachLine;
  		getline(file, eachLine);
		stringstream elementsOfLine(eachLine);
		eachLine.clear();
		vector<double> eachRow;
		while(elementsOfLine)
		{
			elementsOfLine >> data;
			eachRow.push_back(data);
		}	
		if(eachRow.size() > 2) {num_lines++;}
		table_func.push_back(eachRow);
		eachRow.clear();
		elementsOfLine.clear();
	}
	file.close();
	return num_lines;
}
/*************************************************************************************
 ************ 
 ************ Count the number of lines in a file
 ************	
 **************************************************************************************/
int countLinesInFile(const char *fileName)
{
	ifstream file;
	file.open(fileName); 	
	
	int num_lines = 0;
	
	/*read labels*/
	string eachLine;
	getline(file, eachLine);
	eachLine.clear();
	num_lines++;
	
	while(file)
	{
		string eachLine;
  		getline(file, eachLine);
		eachLine.clear();
		num_lines++;
	}
	file.close();
	return num_lines;
}

/*************************************************************************************
 ************ 
 ************	Count number of elements in a given range in radius
 ************	
 **************************************************************************************/
int countElementsInInterval(ifstream &file, int meshElementsRinOut, double rInit, double rFinal, int inwardOutward_func)
{	
	int elementsInRange = 0, i = 0;
	double rFile = 0.0;
	
	file.seekg(0, ios_base::beg);
	while(i <= meshElementsRinOut)
	{
		string eachLine;
		char *stopstring;
		/*read each line into a string*/
  		getline(file, eachLine);
		/*take the string, convert it to double and stop as soon as you hit a white space*/
		rFile = strtod(eachLine.c_str(), &stopstring);
		if(inwardOutward_func == 0 && rFile <= rFinal && rFile >= rInit) //fit to surface
		{
			elementsInRange++;
			if(rFile < rInit){break;}
		}
		if(inwardOutward_func == 1 && rFile >= rFinal && rFile <= rInit) //fit to surface
		{
			elementsInRange++;
			if(rFile > rInit){break;}
		}
		
		eachLine.clear();
		i++;
	}
	file.clear();
	return elementsInRange;
}
/*************************************************************************************
 ************ 
 ************	Fill In the various riccati elements
 ************	
 **************************************************************************************/
void fillInRiccatiElements(ifstream &file, int meshElementsRinOut, double rInit, double rFinal, 
						   int inwardOutward_func, int sizeRij_calc_func, double *csi_calc_func, double **Rij_calc_func, 
						   vector<vector<double> > &table_func)
{
	int i = 0, j = 0, k = 0, elementsInRange = 0;
	double data = 0.0;
	
	file.seekg(0, ios_base::beg);

	while(i <= meshElementsRinOut)
	{
		string eachLine;
  		getline(file, eachLine);
		stringstream elementsOfLine(eachLine);
		eachLine.clear();
		vector<double> eachRow;
		while(elementsOfLine)
		{			
			elementsOfLine >> data;
			eachRow.push_back(data);
		}

		if(inwardOutward_func == 0 && eachRow[0] <= rFinal && eachRow[0] >= rInit) //fit to surface			
		{
			table_func.push_back(eachRow);
			elementsInRange++;
			if(eachRow[0] < rInit){break;}
		}
		if(inwardOutward_func == 1 && eachRow[0] >= rFinal && eachRow[0] <= rInit) //fit to surface
		{
			table_func.push_back(eachRow);
			elementsInRange++;
			if(eachRow[0] > rInit){break;}
		}		
		eachRow.clear();
		elementsOfLine.clear();
		i++;
	}	
	
	/*invert the order of the elements 
	 - from surface to fit becomes from fit to surface
	 - from center to fit becomes from fit to center
	 */
	k = 0;
	for (i = elementsInRange-1; i >= 0 ; i--)
	{
		csi_calc_func[k] = table_func[i][0];
		for (j = 0; j<sizeRij_calc_func; j++) {Rij_calc_func[k][j] = table_func[i][j + 2];}
		k++;
	}
	
	file.clear();
	table_func.clear();	
}
/*************************************************************************************
 ************ 
 ************	Open and label integrator state output files
 ************	
 **************************************************************************************/
void OpenAndLabel_rkf45_state(const char * ActualfileName, int counterPermutation, ofstream &rkf45_state_inOut, int rowsRic_func)
{
	string permutationCntrLabel;        // string which will contain the result					
	ostringstream convert;			    // stream used for the conversion					
	
	convert << counterPermutation;      // insert the textual representation of 'Number' in the characters in the stream					
	permutationCntrLabel = convert.str(); // set 'Result' to the contents of the stream

	string file_rkf45_state	= ActualfileName + permutationCntrLabel + ".dat";						

	rkf45_state_inOut.open(file_rkf45_state.c_str());
	rkf45_state_inOut<< setiosflags(ios_base::scientific);			
	labelF_integratorState_inOut(rkf45_state_inOut,rowsRic_func);	
	
	permutationCntrLabel.clear();
	convert.clear();
}

/*************************************************************************************
 ************ 
 ************	Open and label integrator state output files
 ************	
 **************************************************************************************/
void OpenAndLabel_Rij_inOut(const char * ActualfileName, int counterPermutation, 
							ofstream &elmntsR_inOut, int rowsRic_func, int tidesFlag_func)
{
	string permutationCntrLabel;        // string which will contain the result					
	ostringstream convert;			    // stream used for the conversion					
	
	convert << counterPermutation;      // insert the textual representation of 'Number' in the characters in the stream					
	permutationCntrLabel = convert.str(); // set 'Result' to the contents of the stream
	
	
	string fileRiccati	= ActualfileName + permutationCntrLabel + ".dat";						
	
	elmntsR_inOut.open(fileRiccati.c_str());
	elmntsR_inOut<< setiosflags(ios_base::scientific);			

	labelF_Riccati_inOut(elmntsR_inOut, rowsRic_func, tidesFlag_func);		

	permutationCntrLabel.clear();
	convert.clear();
}

/*************************************************************************************
 ************ 
 ************	Fill In the integrator state form rkf45
 ************	
 **************************************************************************************/
void fillIn_rkf45_state(ifstream &file, int meshElementsRinOut, double rInit, double rFinal, 
						   int inwardOutward_func, int size_rkf45_state_func, double **rRel_and_rkf45_state_func, 
						   vector<vector<double> > &table_func)
{
	int i = 0, j = 0, k = 0, elementsInRange = 0;
	double data = 0.0;
	
	file.seekg(0, ios_base::beg);
	
	while(i <= meshElementsRinOut)
	{
		string eachLine;
  		getline(file, eachLine);
		stringstream elementsOfLine(eachLine);
		eachLine.clear();
		vector<double> eachRow;
		while(elementsOfLine)
		{			
			elementsOfLine >> data;
			eachRow.push_back(data);
		}
		
		if(inwardOutward_func == 0 && eachRow[0] <= rFinal && eachRow[0] >= rInit) //fit to surface			
		{
			table_func.push_back(eachRow);
			elementsInRange++;
			if(eachRow[0] < rInit){break;}
		}
		if(inwardOutward_func == 1 && eachRow[0] >= rFinal && eachRow[0] <= rInit) //fit to surface
		{
			table_func.push_back(eachRow);
			elementsInRange++;
			if(eachRow[0] > rInit){break;}
		}		
		eachRow.clear();
		elementsOfLine.clear();
		i++;
	}	
	
	/*invert the order of the elements 
	 - from surface to fit becomes from fit to surface
	 - from center to fit becomes from fit to center
	 */
			
	k = 0;
	for (i = elementsInRange-1; i >= 0 ; i--)
	{
		for (j = 0; j < size_rkf45_state_func; j++)
			rRel_and_rkf45_state_func[k][j] = table_func[i][j];
		k++;
	}

	file.clear();
	table_func.clear();	
}
/*************************************************************************************
 ************ 
 ************	This routine writes the labels in the file containing
 ************	the Riccati elements Rij at the fitting point.
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write labels
 ************	- rowsRic_func = rows of R to count the riccati elements to label
 ************	
 **************************************************************************************/
void labelF_Riccati_xFit(ofstream &fileName, int rowsRic_func)
{
	int i = 0, j = 0;
	
	fileName  <<"..........omega_r..........log10(omega2)_r...............";
	fileName  <<"omega_i..............";
	for (i = 0; i < rowsRic_func; i++)
		for (j = 0; j < rowsRic_func; j++)
		{
			fileName  << "R" << i << j <<"_o....................";
		}
	for (i = 0; i < rowsRic_func; i++)
		for (j = 0; j < rowsRic_func; j++)
			fileName  << "R" << i << j <<"_i....................";
	
	fileName  <<"detRin_minus_Rout_r........detRin_minus_Rout_i"<< endl;
}
/*************************************************************************************
 ************ 
 ************	This routine writes in a file the Riccati elements at the fitting point
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write Rij
 ************	- omega_func_r = real eigenfrequency
 ************	- omega_func_i = imaginary eigenfrequency
 ************	- matR_in_func = matrix R from inward integration (surface-fit) at fitting point.
 ************	- matR_out_func = matrix R from outward integration (center-fit) at fitting point.
 ************	- vecDetDiff_r_i_func = vector containing the riccati condition at the fitting point:
 ************	(element#0: det(R_in - R_out) real,  element#1: det(R_in - R_out) imaginary) 
 ************	- rowsRic_func = rows of matrix R
 ************	- sizeDetVector_func = size of the vector vecDetDiff_r_i_func
 ************	
 **************************************************************************************/
void writeF_Riccati_xFit(ofstream &fileName, double omega_func_r, double omega_func_i, gsl_matrix *matR_in_func, gsl_matrix *matR_out_func, 
						 gsl_vector *vecDetDiff_r_i_func, int rowsRic_func, int sizeDetVector_func)
{
	int i = 0, j = 0;
	
	fileName << setprecision(16) << omega_func_r << "  " 
	<< log10(omega_func_r*omega_func_r) << "   ";
	
	fileName << setprecision(16) << omega_func_i << "  ";
	
	for (i = 0; i< rowsRic_func; i++)
		for (j = 0; j<rowsRic_func; j++)	
			fileName << gsl_matrix_get(matR_out_func,i,j) << "  ";
	
	for (i = 0; i< rowsRic_func; i++)
		for (j = 0; j<rowsRic_func; j++)	
			fileName << gsl_matrix_get(matR_in_func,i,j) << "  ";
	
	for(i = 0; i<sizeDetVector_func; i++)
		fileName << gsl_vector_get(vecDetDiff_r_i_func, i) << "  ";
	
	fileName << endl;
}

/*************************************************************************************
 ************ 
 ************	This routine writes the labels in the files containing
 ************	the Riccati elements calculated during the inward and
 ************	outward integration (center - fit and surface-fit)
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write labels
 ************	- rowsRic_func = rows of R to count the riccati elements to label
 ************	
 **************************************************************************************/
void labelF_Riccati_inOut(ofstream &fileName, int rowsRic_func, int tidesFlag_func)
{
	int i = 0, j = 0;
	
	fileName  <<"...............csiRel...............csi...................";
	for (i = 0; i < rowsRic_func; i++)
		for (j = 0; j < rowsRic_func; j++)
			fileName  << "r" << i << j <<"......................";
	
	fileName  << "omega_re.......................omega_im.......................l";
	if(tidesFlag_func)
	{
		fileName  << "Pspin(min).......................Porb(min).......................";
		fileName  << "m.......................k";
	}

	fileName  << endl;

}
/*************************************************************************************
 ************ 
 ************	This routine writes in a file the Riccati elements calculated during 
 ************	the inward and outward integration (center - fit and surface-fit)
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write Rij
 ************	- csiRel_func = mesh point r dimensionless
 ************	- csiAbs_func = mesh point r with dimensions
 ************	- vec_y_Riccati_all_func = vector containing the Rij to be written
 ************	- sizeR_func = size of the riccati matrix to count the number of elements to write
 ************	- omega2re_func = real eigenfrequency
 ************	- omega2im_func = imaginary eigenfrequency
 ************	- deg_l_func = l-order of the mode.
 ************	
 **************************************************************************************/
void writeF_Riccati_inOut(ofstream &fileName, double csiRel_func, double csiAbs_func, 
						 gsl_vector *vec_y_Riccati_all_func, 
						 int sizeR_func, double omega_r_func, double omega_i_func, double deg_l_func, 						  
						  double deg_m_func, double deg_k_func, double Pspin_minutes_func, double Porb_minutes_func, 
						  int tidesFlag_func)
{
	int i = 0;
	fileName << setprecision(16) << csiRel_func << "  " << csiAbs_func << "   "; 

	for (i = 0; i<sizeR_func; i++)
		fileName << gsl_vector_get(vec_y_Riccati_all_func,i) << "   " ;

	fileName << omega_r_func << "   " << omega_i_func << "   " <<  deg_l_func ;
	
	if(tidesFlag_func)
	{
		fileName << "   " << Pspin_minutes_func << "   " << Porb_minutes_func << "   ";		
		fileName << deg_m_func << "   " << deg_k_func;		
	}
	
	fileName << endl;

}
/*************************************************************************************
 ************ 
 ************	This routine writes the labels in the files containing
 ************	the eigenfunctions.
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write labels
 ************	- rowsRic_func = rows of R to count the number of eigenfunctions
 ************	- adiabatic_func 
 ************	- tidesFlag_func = flag for tides if I need y7 ad y8
 ************	
 **************************************************************************************/

void labelF_eigenfunctions(ofstream &fileName, int rowsRic_func, int adiabatic_func, int tidesFlag_func)
{
	
	int i = 0, j = 0, sizeUplusV = rowsRic_func+rowsRic_func;
	
	fileName  <<"........csiRel....................csi......................";

	if(adiabatic_func)		
	{
		if(tidesFlag_func)
			for (i = 1; i <= sizeUplusV; i++)			
			{	
				j = i;
 				if(i == 5 || i == 6){j = i+2;}
					fileName  << "y" << j <<".....................";				
			}
		else
			for (i = 1; i <= sizeUplusV; i++)
				fileName  << "y" << i <<".....................";

	}//if(adiabatic_func)		
	else
	{
		for (i = 1; i <= rowsRic_func; i++)
			fileName  << "y" << i <<"real..................";
		for (i = 1; i <= rowsRic_func; i++)
			fileName  << "y" << i <<"imag..................";
	}
		
	for (i = 0; i< rowsRic_func; i++)
		for (j = 0; j<rowsRic_func; j++)
			fileName  << "r'" <<i<<j<<"......................";

	fileName  <<"....Omega_r........................Omega_i......................l";
	fileName  << ".....................rho..........................g..........";

	if(tidesFlag_func)
	{
		fileName  << "...........Pspin(min).............Porb(min)";
		fileName  << "..................m..................k.........";
		fileName  << "a/dadt_tides (yr)..........a/dadt_GR (yr)..........omega/dOmegadt_GR (yr)";
		fileName  << "..........argFlmk(deg)..........modFlmk";

	}
	
	fileName<<endl;

	
}
/*************************************************************************************
 ************ 
 ************	This routine writes in a file the eigenfunctions
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write Rij
 ************	- csiRel_func = mesh point r dimensionless
 ************	- csiAbs_func = mesh point r with dimensions
 ************	- omega2_r_func = real eigenfrequency
 ************	- omega2_i_func = imaginary eigenfrequency
 ************	- deg_l_func = l-order of the mode.
 ************	- rho_func = density at a given mesh point
 ************	- g_func = gravity at a given mesh point
 ************	- vecUoriginal_func = Riccati vector containing part of the eigenfunctions (see below)
 ************	- vecVoriginal_func = Riccati vector containing part of the eigenfunctions (see below)
 ************	- matRprime_func = matrix containing the permuted Riccati components
 ************	- rowsRic_func = rows of the riccati matrix to count the number of elements to write
 ************	
 ************	
 ************	
 ************	Adiabatic + no tides case:
 ************	u0 = y1, u1 = y2, v0 = y3, v1 = y4 
 ************	
 ************	Adiabatic + tides case:
 ************	u0 = y1, u1 = y2, u2 = y7
 ************	v0 = y3, v1 = y4, v2 = y8 
 ************	
 ************	Non adiabatic case:
 ************	u0 = y1r, u1 = y2r, u2 = y5r, u3 = y1i, u4 = y2i, u5 = y5i, 
 ************	v0 = y3r, v1 = y4r, v2 = y6r, v3 = y3i, v4 = y4i, v5 = y6i 
 ************	
 **************************************************************************************/
int writeF_eigenfunctions(ofstream &fileName, double csiRel_func, double csiAbs_func, 
						  double omega_r_func, double omega_i_func, double deg_l_func, double rho_func, double g_func,
						  gsl_vector *vecUoriginal_func, gsl_vector *vecVoriginal_func, 
						  gsl_matrix *matRprime_func, int rowsRic_func, int adiabatic_func, int tidesFlag_func, 
						  double Pspin_func, double Porb_func, double deg_m_func ,double deg_k_func, 
						  double tidesTimescale_a_func, double GRtimescale_a_func, double tidesTimescale_Omega_func, 
						  double arg_Flmk_func, double Mod_Flmk_func)
{

	int i = 0, j = 0, columnsEigenfuncFile = 0;
	
	fileName << setprecision(16) << csiRel_func << "  ";
	columnsEigenfuncFile++;
	fileName << csiAbs_func << "   ";
	columnsEigenfuncFile++;
	
	if (adiabatic_func)
	{
		if(tidesFlag_func)
		{
			fileName << gsl_vector_get(vecUoriginal_func,0) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecUoriginal_func,1) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,0) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,1) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecUoriginal_func,2) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,2) << "  ";			
			columnsEigenfuncFile++;
		}//if(tidesFlag_func)
		else
		{
			for (i = 0; i < rowsRic_func; i++)
			{
				fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			for (i = 0; i < rowsRic_func; i++)
			{
				fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;			
				columnsEigenfuncFile++;
			}
		}//else
		for (i = 0; i < rowsRic_func; i++)
			for (j = 0; j < rowsRic_func; j++)
			{
				fileName << gsl_matrix_get(matRprime_func,i,j) << "  " ;	
				columnsEigenfuncFile++;
			}
		omega_i_func = 0;
	}//if(adiabatic)
	else
	{
		
		for (i = 0; i <=1; i++)
		{
			fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
			columnsEigenfuncFile++;
		}
		
		for (i = 0; i <=1; i++)
		{
			fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;
			columnsEigenfuncFile++;
		}
		
		
		fileName << gsl_vector_get(vecUoriginal_func,2) << "  " ;
		columnsEigenfuncFile++;
		fileName << gsl_vector_get(vecVoriginal_func,2) << "  " ;
		columnsEigenfuncFile++;

		if(tidesFlag_func)
		{
			fileName << gsl_vector_get(vecUoriginal_func,3) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,3) << "  " ;
			columnsEigenfuncFile++;

			for (i = 4; i <=5; i++)
			{
				fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			for (i = 4; i <=5; i++)
			{
				fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			fileName << gsl_vector_get(vecUoriginal_func,6) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,6) << "  " ;
			columnsEigenfuncFile++;

			fileName << gsl_vector_get(vecUoriginal_func,7) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,7) << "  " ;
			columnsEigenfuncFile++;
		}//if(tidesFlag_func)
		else
		{
			for (i = 3; i <=4; i++)
			{
				fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			for (i = 3; i <=4; i++)
			{
				fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}			
			fileName << gsl_vector_get(vecUoriginal_func,5) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,5) << "  " ;
			columnsEigenfuncFile++;
		
		}//else no tides
		
		for (i = 0; i < rowsRic_func; i++)
			for (j = 0; j < rowsRic_func; j++)
			{
				fileName << gsl_matrix_get(matRprime_func,i,j) << "  " ;		
				columnsEigenfuncFile++;
			}
		
	}

	fileName << omega_r_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << omega_i_func <<  "   " ;
	columnsEigenfuncFile++;
	fileName << deg_l_func << "   ";
	columnsEigenfuncFile++;
	fileName << rho_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << g_func;	
	columnsEigenfuncFile++;

	if(tidesFlag_func)
	{
		fileName << "   " << Pspin_func <<   "   ";
		columnsEigenfuncFile++;
		fileName << Porb_func << "   ";
		columnsEigenfuncFile++;
		fileName << deg_m_func << "   ";
		columnsEigenfuncFile++;
		fileName << deg_k_func << "   ";
		columnsEigenfuncFile++;
		fileName << tidesTimescale_a_func << "   ";
		columnsEigenfuncFile++;
		fileName << GRtimescale_a_func << "   ";
		columnsEigenfuncFile++;
		fileName << tidesTimescale_Omega_func << "   ";
		columnsEigenfuncFile++;
		fileName << arg_Flmk_func << "   ";
		columnsEigenfuncFile++;
		fileName << Mod_Flmk_func;	
		columnsEigenfuncFile++;
	}

	fileName << endl;
	return columnsEigenfuncFile;
}

/*************************************************************************************
 ************ 
 ************	This routine writes the labels in the files containing
 ************	the eigenfunctions without the riccati coefficients.
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write labels
 ************	- rowsRic_func = rows of R to count the number of eigenfunctions
 ************	- adiabatic_func 
 ************	- tidesFlag_func = flag for tides if I need y7 ad y8
 ************	
 **************************************************************************************/
void labelF_eigenfunctions_and_unperturbed_model(ofstream &fileName, int rowsRic_func, int adiabatic_func, int tidesFlag_func)
{
	
	int i = 0, j = 0, sizeUplusV = rowsRic_func+rowsRic_func;
	
	fileName  <<"All quantities are dimensionless as given by the stellar model (G = M = R = Rgas = 1), unless specified" << endl;

	fileName  <<"........csiRel.............";
	fileName  <<"csiAbs (Rsun)...................";
	
	if(adiabatic_func)		
	{
		if(tidesFlag_func)
			for (i = 1; i <= sizeUplusV; i++)			
			{	
				j = i;
 				if(i == 5 || i == 6){j = i+2;}
				fileName  << "y" << j <<"......................";				
			}
		else
			for (i = 1; i <= sizeUplusV; i++)
				fileName  << "y" << i <<"......................";
		
	}//if(adiabatic_func)		
	else
	{
		for (i = 1; i <= rowsRic_func; i++)
			fileName  << "y" << i <<"real..................";
		for (i = 1; i <= rowsRic_func; i++)
			fileName  << "y" << i <<"imag..................";
	}

	fileName  << "Pspin(min)................";
	fileName  << "rho......................";
	fileName  << "g..........................";
	fileName  << "P......................";
	fileName  << "Gamma1................";
	fileName  << "N2.....................";
	fileName  << "delAd......................";
	fileName  << "T (g/mol)................";
	fileName  << "cP(mol/g)....................";
	fileName  << "entropy................";
	fileName  << "Lr.......................";
	
	if(tidesFlag_func)
	{
		fileName  << "a/dadt_T (yr)..........e/dedt_T (yr)..........omega/dOmegadt_T(yr)";
		fileName  << "..........a/dadt_GR (yr)..........argFlmk(deg)..........modFlmk";		
		fileName  << "..........linearityViolation (1/0)";		
		fileName  << "...........(1/J)dJ/dt";		
		fileName  << ".............JdotSpin/Jspin";		
		fileName  << ".............JdotOrb/Jorb";		
	}
	fileName  << ".................V..........";

	fileName<<endl;
}


/*******************************************************************************************
 ************ 
 ************	This routine writes in a file the eigenfunctions without the Riccati coeffs.
 ************	and with some of the model's unperturbed properties
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write Rij
 ************	- vecUoriginal_func = Riccati vector containing part of the eigenfunctions (see below)
 ************	- vecVoriginal_func = Riccati vector containing part of the eigenfunctions (see below)
 ************	- Desired quantities
 ************	
 ************	Adiabatic + no tides case:
 ************	u0 = y1, u1 = y2, v0 = y3, v1 = y4 
 ************	
 ************	Adiabatic + tides case:
 ************	u0 = y1, u1 = y2, u2 = y7
 ************	v0 = y3, v1 = y4, v2 = y8 
 ************	
 ************	Non adiabatic case:
 ************	u0 = y1r, u1 = y2r, u2 = y5r, u3 = y1i, u4 = y2i, u5 = y5i, 
 ************	v0 = y3r, v1 = y4r, v2 = y6r, v3 = y3i, v4 = y4i, v5 = y6i 
 ************	
 **************************************************************************************/
int writeF_eigenfunctions_and_unperturbed_model(ofstream &fileName, double csiRel_func, double csiAbs_func, 						  
												gsl_vector *vecUoriginal_func, gsl_vector *vecVoriginal_func, 
												double rho_func, double g_func, double P_func, double Gamma1_func, 
												double N2_func, double delAd_func, double T_func, double Cp_func, 
												double S_func, double Lr_func, double tidesTimescale_a_func, double tidesTimescale_e_func, 
												double tidesTimescale_Omega_func, double GRtimescale_a_func, 
												double arg_Flmk_func, double Mod_Flmk_func, double Pspin_func, int rowsRic_func, 
												int adiabatic_func, int tidesFlag_func, int linearityViol_func, 
												double Jorb_func, double Jspin_func, double JdotOrb_func, double JdotSpin_func, double V_func)
{
	
	int i = 0, columnsEigenfuncFile = 0;
	
	fileName << setprecision(16) << csiRel_func << "  ";
	columnsEigenfuncFile++;
		fileName << csiAbs_func << "   ";
	columnsEigenfuncFile++;
	
	if (adiabatic_func)
	{
		if(tidesFlag_func)
		{
			fileName << gsl_vector_get(vecUoriginal_func,0) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecUoriginal_func,1) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,0) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,1) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecUoriginal_func,2) << "  ";
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,2) << "  ";			
			columnsEigenfuncFile++;
		}//if(tidesFlag_func)
		else
		{
			for (i = 0; i < rowsRic_func; i++)
			{
				fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			for (i = 0; i < rowsRic_func; i++)
			{
				fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;			
				columnsEigenfuncFile++;
			}
		}//else
	}//if(adiabatic)
	else
	{
		
		for (i = 0; i <=1; i++)
		{
			fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
			columnsEigenfuncFile++;
		}
		
		for (i = 0; i <=1; i++)
		{
			fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;
			columnsEigenfuncFile++;
		}
		
		
		fileName << gsl_vector_get(vecUoriginal_func,2) << "  " ;
		columnsEigenfuncFile++;
		fileName << gsl_vector_get(vecVoriginal_func,2) << "  " ;
		columnsEigenfuncFile++;
		
		if(tidesFlag_func)
		{
			fileName << gsl_vector_get(vecUoriginal_func,3) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,3) << "  " ;
			columnsEigenfuncFile++;
			
			for (i = 4; i <=5; i++)
			{
				fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			for (i = 4; i <=5; i++)
			{
				fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			fileName << gsl_vector_get(vecUoriginal_func,6) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,6) << "  " ;
			columnsEigenfuncFile++;
			
			fileName << gsl_vector_get(vecUoriginal_func,7) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,7) << "  " ;
			columnsEigenfuncFile++;
		}//if(tidesFlag_func)
		else
		{
			for (i = 3; i <=4; i++)
			{
				fileName << gsl_vector_get(vecUoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}
			
			for (i = 3; i <=4; i++)
			{
				fileName << gsl_vector_get(vecVoriginal_func,i) << "  " ;
				columnsEigenfuncFile++;
			}			
			fileName << gsl_vector_get(vecUoriginal_func,5) << "  " ;
			columnsEigenfuncFile++;
			fileName << gsl_vector_get(vecVoriginal_func,5) << "  " ;
			columnsEigenfuncFile++;
			
		}//else no tides		
	}
	fileName << Pspin_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << rho_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << g_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << P_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << Gamma1_func <<  "   ";	
	columnsEigenfuncFile++;
	fileName << N2_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << delAd_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << T_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << Cp_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << S_func <<  "   ";
	columnsEigenfuncFile++;
	fileName << Lr_func << "   ";
	columnsEigenfuncFile++;

	if(tidesFlag_func)
	{		
		fileName << "   " << tidesTimescale_a_func << "   ";
		columnsEigenfuncFile++;
		fileName << tidesTimescale_e_func << "   ";
		columnsEigenfuncFile++;
		fileName << tidesTimescale_Omega_func << "   ";
		columnsEigenfuncFile++;
		fileName << GRtimescale_a_func << "   ";
		columnsEigenfuncFile++;
		fileName << arg_Flmk_func << "   ";
		columnsEigenfuncFile++;
		fileName << Mod_Flmk_func << "          ";
		columnsEigenfuncFile++;
		fileName << linearityViol_func << "          ";
		columnsEigenfuncFile++;
		fileName << (JdotOrb_func + JdotSpin_func)/(Jorb_func + Jspin_func) << "   ";
		columnsEigenfuncFile++;
		fileName << JdotSpin_func/Jspin_func << "   ";
		columnsEigenfuncFile++;
		fileName << JdotOrb_func/Jorb_func << "   ";
		columnsEigenfuncFile++;
	}
	fileName << V_func;	
	columnsEigenfuncFile++;

	fileName << endl;
	return columnsEigenfuncFile;
}
/*************************************************************************************
 ************ 
 ************	This routine writes the labels in the files containing
 ************	the global properties of the system
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write labels
 ************	
 **************************************************************************************/
void labelF_globalProperties(ofstream &fileName)
{
	
	fileName  <<".........M1(Msun)";
	fileName  <<"...............R1(Rsun)";
	fileName  <<"...............L1(Lsun)";
	fileName  <<"..................Teff1(K)";
	fileName  <<"...............I1(MsunRsun^2)";
	fileName  <<"...............M2(Msun)";
	fileName  <<"...............Pspin(min)";
	fileName  <<"...............Porb(min)";
	fileName  <<"...............a(Rsun)";
	fileName  <<"...............ecc";
	fileName  <<"........................omega_R";
	fileName  <<".................omega_I";
	fileName  <<"...............PbreakUp(min)";
	fileName  <<"....................l";
	fileName  <<".......................m";
	fileName  <<".......................k";
	fileName  <<".......................epsilonT";
	fileName  <<"....................clmk";
	fileName  <<"..................G2lmk";
	fileName  <<"..................G3lmk";
	fileName  <<"..................klmk";
	
	fileName  << endl;
	
}

/*************************************************************************************
 ************ 
 ************	This routine writes in a file the global properties of the star.
 ************	
 ************	Input:
 ************	- fileName = ofstream file where to write Rij
 ************	- global properties to be written
 ************	
 **************************************************************************************/
void writeF_globalProperties(ofstream &fileName, double M1_Msun_func, double R1_Rsun_func,  
							 double L1_Lsun_func, double Teff_K_func, double momInertia_MsunRsun2_func,  
							 double M2_Msun_func, double Pspin_minutes_func, double Porb_minutes_func,  
							 double a_Rsun_func, double ecc_func, double omega_r_func, double omega_i_func, 
							 double PbreakUp_min_func, double l_func, double deg_m_func, double deg_k_func, 
							 double epsilon_T_func, double Clmk_func, double Glmk_2_func, double Glmk_3_func, 
							 double klmk_func)
{
	
	fileName << setprecision(16) << M1_Msun_func << "   " << R1_Rsun_func << "   " << L1_Lsun_func << "   " 
	<< Teff_K_func << "   " << momInertia_MsunRsun2_func << "   " << M2_Msun_func << "   " << Pspin_minutes_func << "   " 
	<< Porb_minutes_func << "   " << a_Rsun_func << "   " << ecc_func << "   " << omega_r_func << "   " << omega_i_func << "   " 
	<< PbreakUp_min_func << "   " << l_func << "   " << deg_m_func << "   " << deg_k_func << "   " << epsilon_T_func << "   " 
	<< Clmk_func << "   " << Glmk_2_func << "   " << Glmk_3_func << "   " << klmk_func << endl;
}
/*************************************************************************************
 ************ 
 ************	This routine writes the labels in the files containing
 ************	rkf45 state. rkf45 state can be used to compute Rij during the 
 ************	calculation of the eigenfunctions
 ************ 
 **************************************************************************************/
void labelF_integratorState_inOut(ofstream &fileName, int rowsRic_func)
{
	int i = 0, j = 0;
	
	fileName  <<"...........csiRel..............";
	for (i = 0; i < rowsRic_func; i++)
		for (j = 0; j < rowsRic_func; j++)
			fileName  << "k1_r" << i << j <<"...................k6_r" << i << j <<"....................";
	
	fileName  << endl;
	
}
/*************************************************************************************
 ************	This routine writes the files containing rkf45 state
 **************************************************************************************/
void writeF_integratorState_inOut(ofstream &fileName, gsl_vector *vecIntegratorState_func, int size_rkf45_state_func)
{
	for (int i = 0; i < size_rkf45_state_func; i++)
		fileName << setprecision(16) << gsl_vector_get(vecIntegratorState_func, i) << "   " ;
	
	fileName << endl;
}

/*************************************************************************************
 ************ 
 ************	This routine writes an error message and exits CAFein
 ************	
 **************************************************************************************/

void errorMessageAndExit(const char * fileName, const char * errorMessage)
{

	cerr<< "#####################################################################"<< endl;;
	cerr<< "######## "<< endl;
	cerr<< "######## ERROR IN: "<<  fileName     << endl;   
	cerr<< "######## " << errorMessage  << endl;
	cerr<< "######## "<< endl;
	cerr<< "######################################################################"<< endl;;
	exit (1);			
	
}








