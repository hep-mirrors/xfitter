#include "LUPinverse.h"
#include <math.h>

int LUPDecompose(cl2DArray<double> *A, int N, double Tol, int *Pdec)
{
  int i, j, k, imax; 
  double maxA, ptr, absA;
   
  for (i = 0; i <= N; i++)
    Pdec[i] = i; // Unit permutation matrix, Pdec[N] initialized with N

  for (i = 0; i < N; i++) {
    maxA = 0.0;
    imax = i;

    for (k = i; k < N; k++)
      if ((absA = fabs((*A)(k,i))) > maxA) { 
        maxA = absA;
        imax = k;
      }

    if (maxA < Tol) return 1; // failure, matrix is degenerate

    if (imax != i) {
      // pivoting Pdec
      j = Pdec[i];
      Pdec[i] = Pdec[imax];
      Pdec[imax] = j;

      // pivoting rows of A
      for (int m = 0; m < N; m++) {
        ptr = (*A)(i,m);
        (*A)(i,m) = (*A)(imax,m);
        (*A)(imax,m) = ptr;
      }

      // counting pivots starting from N (for determinant)
      Pdec[N]++;
    }

    for (j = i + 1; j < N; j++) {
      (*A)(j,i) /= (*A)(i,i);
      for (k = i + 1; k < N; k++)
        (*A)(j,k) -= (*A)(j,i) * (*A)(i,k);
    }
  }

  return 0; // decomposition done 
}

void LUPSolve(cl2DArray<double> *A, int *Pdec, double *b, int N, double *x)
{
  for (int i = 0; i < N; i++) {
    x[i] = b[Pdec[i]];
    for (int k = 0; k < i; k++)
      x[i] -= (*A)(i,k) * x[k];
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int k = i + 1; k < N; k++)
      x[i] -= (*A)(i,k) * x[k];
    x[i] /= (*A)(i,i);
  }
}

void LUPInvert(cl2DArray<double> *A, int *Pdec, int N, cl2DArray<double> *IA)
{
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      (*IA)(i,j) = (Pdec[i] == j) ? 1.0 : 0.0;
      for (int k = 0; k < i; k++)
        (*IA)(i,j) -= (*A)(i,k) * (*IA)(k,j);
    }

    for (int i = N - 1; i >= 0; i--) {
      for (int k = i + 1; k < N; k++)
        (*IA)(i,j) -= (*A)(i,k) * (*IA)(k,j);
      (*IA)(i,j) /= (*A)(i,i);
    }
  }
}

double LUPDeterminant(cl2DArray<double> *A, int *Pdec, int N)
{
  double det = (*A)(0,0);
  for (int i = 1; i < N; i++)
    det *= (*A)(i,i);
  return (Pdec[N] - N) % 2 == 0 ? det : -det;
}
