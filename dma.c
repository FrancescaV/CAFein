#include<stdio.h>
#include<stdlib.h>
#include <iostream>
#include "dma.h"
using namespace std;

/*************************************************************************************
************** This routine is to dynamically allocate memory:
************** For every dimension, the arrays are made one entry bigger to allow 
************** for an empty, unused slot with index zero 
**************************************************************************************/


float *ffun(int n){
	n=n+1;
	float *temp;
	temp = (float*) malloc (n* sizeof(float));
	if (temp==NULL) {cout << "ERROR-----Memory For (float *) Not Available" << endl;}
	return temp;
}

double *dfun(int n){
	n=n+1;
	double *tempo;
	tempo =(double*) malloc(n* sizeof(double));
	if (tempo==NULL) {cout << "ERROR-----Memory For (double *) Not Available" << endl;}
	return tempo;
}

int *ifun(int n){
	n=n+1;
	int *tempp;
	tempp = (int*) malloc(n* sizeof(int));
	if (tempp==NULL) {cout << "ERROR-----Memory For (int *) Not Available" << endl;}
	return tempp;
}


float **ffunc(int n, int m){
	int i;
	float **temp;
	n=n+1;
	m=m+1;
	temp=(float**) malloc (n* sizeof(float));
	for(i=0; i<n; i++){
		temp[i]=(float*) malloc (m* sizeof(float));
	}
	if (temp==NULL) {cout << "ERROR-----Memory For (float **) Not Available" << endl;}
	return temp;
}

double **dfunc(int n, int m){
	int i;
	double **tempo;
	n=n+1;
	m=m+1;
	tempo=(double**) malloc (n* sizeof(double));
	for(i=0; i<n; i++){
		tempo[i]=(double*) malloc (m* sizeof(double));
	}
	if (tempo==NULL) {cout << "ERROR-----Memory For (double **) Not Available" << endl;}
	return tempo;
}

  
int **ifunc(int n, int m){
	int i;
	int **tempp;
	n=n+1;
    m=m+1;
	tempp=(int**) malloc (n* sizeof(int));
	for(i=0; i<n; i++){
		tempp[i]=(int*) malloc (m* sizeof(int));
	}
	if (tempp==NULL) {cout << "ERROR-----Memory For (int **) Not Available" << endl;}

	return tempp;
}

