 
/*
   @file ReactiontestZMVFNS.cc
   @date 2016-12-08
   @author  AddReaction.py
   Created by  AddReaction.py on 2016-12-08
*/

#include "ReactiontestZMVFNS.h"
#include "iostream"

// the class factories
extern "C" ReactiontestZMVFNS* create() {
  return new ReactiontestZMVFNS();
}

// define QCDNUM function:
extern "C" {
  void zmstfun_(int *id, double *key, double *x, double *q2, double *sf, int *np, int *flag);
}



// Initialize at the start of the computation
void ReactiontestZMVFNS::initAtStart(const string &s)
{
}

// Main function to compute results at an iteration
int ReactiontestZMVFNS::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  // Example implementation of low Q2 reduced cross section:

  // Get bin arrays, check that Q2, x and y are present:

  std::valarray<double> *q2p  = GetBinValues(1,"Q2");  
  std::valarray<double> *xp  = GetBinValues(1,"x");  
  std::valarray<double> *yp  = GetBinValues(1,"y");  
  
  if (q2p == NULL || xp == NULL || yp == NULL ) {
    std::cout << "\n\nFATAL ERROR: DIS NC requires x,Q2 and y bins to be present !!!" << std::endl;
    std::cout << "CHECK THE DATAFILE !!!" << std::endl;
    return 1;
  }
  std::valarray<double> q2 = *q2p;
  std::valarray<double> x = *xp;
  std::valarray<double> y = *yp;


  // Compute u and d-type F2s:

  int Npnt = q2.size();
  std::valarray<double> f2u(Npnt), f2d(Npnt), f2(Npnt);
  std::valarray<double> flu(Npnt), fld(Npnt), fl(Npnt);
  int flag = 0;

  double CNEP2F[] = {0.,0.,1.,0.,1.,0.,0.,0.,1.,0.,1.,0.,0.}; //u
  double CNEM2F[] = {0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.,1.,0.}; //d

  // u-type:
  int id   = 1;
  zmstfun_(&id,&CNEP2F[0], &x[0], &q2[0], &flu[0], &Npnt, &flag);
  id = 2;
  zmstfun_(&id,&CNEP2F[0], &x[0], &q2[0], &f2u[0], &Npnt, &flag);

  // d-type:
  id   = 1;
  zmstfun_(&id,&CNEM2F[0], &x[0], &q2[0], &fld[0], &Npnt, &flag);
  id = 2;
  zmstfun_(&id,&CNEM2F[0], &x[0], &q2[0], &f2d[0], &Npnt, &flag);

  // valarrays are smart enough for that:
  f2 = 4./9.*f2u + 1./9.*f2d;
  fl = 4./9.*flu + 1./9.*fld;

  
  // compute reduced x-section:
  for (int i = 0; i<Npnt; i++) {
    double yplus = 1+(1-y[i])*(1-y[i]);
    double y2    = y[i]*y[i];
    val[i] = f2[i] - y2/yplus * fl[i];
  }  
  return 0;
}

