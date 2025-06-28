/*
Description: Functions for LU decomposition of square matrices
*/

#ifndef LUPINVERSE_H
#define LUPINVERSE_H

#include <iostream>
#include "cl2DArray.h"

// LU Decomposition function declarations
int LUPDecompose(cl2DArray<double> *A, int N, double Tol, int *Pdec);
void LUPSolve(cl2DArray<double> *A, int *Pdec, double *b, int N, double *x);
void LUPInvert(cl2DArray<double> *A, int *Pdec, int N, cl2DArray<double> *IA);
double LUPDeterminant(cl2DArray<double> *A, int *Pdec, int N);

#endif // LUPINVERSE_H
