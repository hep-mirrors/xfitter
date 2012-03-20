/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2005--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#include <cstdlib>
#include <cmath>
#include <cstring>
/*
#include <stdio.h>
#include <float.h>
#include <iostream.h>
*/

using namespace std;

#include "emsg.h"
#include "tblipol.h"
#include "MDipol.h"

//*************************************************************************
static void MDlocate(real_type *x, real_type *xx[], int *npt, int dim, int *ix) {
  int dd;
  for(dd=0; dd < dim; dd++) {
    require(npt[dd] > 1, "MDlocate: too little points to interpolate.");
    ix[dd] = FindIndexR(xx[dd], npt[dd], x[dd]);
    require(ix[dd] >= 0 && ix[dd] < (npt[dd]-1),
      "MDlocate: x[%d] = %g out of range.", dd, x[dd]);
  }
}

//*************************************************************************
static int Index(const int *ind, int *npt, int dim){
  if(dim==1) return ind[0];
  return ind[dim-1] + npt[dim-1]*Index(ind, npt, dim-1);
}

#ifdef TEST
  cout << exp(x[0]) << ", " << exp(x[1]) << endl;
  //cout << npt[0] << ", " << npt[1] << endl;
  cout << ix[0] << ", " << ix[1] << endl;
  cout << exp(xx[1][ix[1]]) << ", " << exp(xx[1][ix[1]+1]) << endl;
  cout << "y0 = " << y0 << "\n----------------" << endl;
  //cin.ignore(100,'\n');
#endif

//*************************************************************************
static real_type Psi(int n, int dim, int *npt, const int *ix, real_type lam[], real_type yy[]){
  if(n==0) return yy[Index(ix,npt,dim)];
  //int *k0 = new int[dim];
  int *ix1 = new int[dim];
  //memcpy(k0,k,dim*sizeof(int));
  memcpy(ix1,ix,dim*sizeof(int));
  //k0[n] = 0;
  ix1[n-1]++;
  real_type v =
    lam[n-1]*Psi(n-1,dim,npt,ix1,lam,yy) + (1-lam[n-1])*Psi(n-1,dim,npt,ix,lam,yy);
  delete[] ix1;
  return v;
}

//*************************************************************************
real_type LinIpol(real_type *x, real_type *xx[], real_type yy[],
          int *npt, int dim/* =1*/) {
//--- linear interpolation in arbitrary # of dimensions
// yy[n1][n2]...[n_dim]
// i[dim] + n[dim]*(i[dim-1] + n[dim-1]*i[dim-2] + ...)
real_type val;
int *ix = new int[dim];
  MDlocate(x, xx, npt, dim, ix);
  //---  x --> lambda
  real_type* lam = new real_type[dim];
  for(int n=0; n < dim; n++)
    lam[n] = (x[n]-xx[n][ix[n]])/(xx[n][ix[n]+1]-xx[n][ix[n]]);
  val = Psi(dim,dim,npt,ix,lam,yy);
  delete[] ix;
  delete[] lam;
  return val;
}
