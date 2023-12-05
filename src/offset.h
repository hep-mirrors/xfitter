/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2013
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _OFFSET_H_
#define _OFFSET_H_

// ------------------------
// #define DEBUG_OFFSET
// ------------------------

#ifdef DEBUG_OFFSET
  #define TNT_BOUNDS_CHECK
#endif

#include <cstdio>
#include <fstream>
#include <cstring>
// #include "matrix.h"
#include "CovMatrix.h"

// ooooooooooooooooooooooooooo
class OffsetCalc_t {
public:
  MATRIX fitparM, fitparP;
  VECTOR fitpar0;
  const int nCosy;
  const int nVarPar;
  
  // ============================================================
  OffsetCalc_t(int n_Cosy, int n_VarPar) : nCosy(n_Cosy), nVarPar(n_VarPar) {
    // --- create matrix for shifted parameters
    fitparM.newsize(nCosy, nVarPar);
    fitparP.newsize(nCosy, nVarPar);
    fitpar0.newsize(nVarPar);
    #ifdef DEBUG_OFFSET
      cout << "--> OffsetCalc_t: " << "nCosy = "<< nCosy <<", nVarPar = "<< nVarPar << endl;
    #endif
  }
  
  // ============================================================
  void NewParams(const double* p, int mu) {
    #ifdef DEBUG_OFFSET
      cout << "--> NewParams, mu = " << mu << endl;
      for(int j=0; j < nVarPar; j++) cout << j <<": "<< p[j] << endl;
    #endif
    int amu = abs(mu);
    if(mu > 0) memcpy(fitparP[amu-1], p, nVarPar*sizeof(double));
    else if(mu < 0) memcpy(fitparM[amu-1], p, nVarPar*sizeof(double));
    else memcpy(fitpar0, p, nVarPar*sizeof(double));
  }
  
  // ============================================================
  MATRIX GetCov() {
    MATRIX DA = mult(fitparP - fitparM, 0.5);
    #ifdef DEBUG_OFFSET
      cout << "--> GetCov, DA = " << DA << endl;
      cout << "      " << transpose_mult(DA,DA) << endl;
    #endif
    return transpose_mult(DA,DA);
  }
  
  // ============================================================
  int FillFromFiles(const string& outdir="") {
    string fpath(outdir);
    if(!fpath.empty() && *(fpath.end()-1) != '/') fpath.append("/");
    ifstream dat;

    fpath.append("params_");
    char ts[16];
    VECTOR v;
    for(int mu = -nCosy; mu <= nCosy; mu++) {
      if(mu) {
	snprintf(ts, sizeof(ts), "%03d%c.txt", abs(mu), mu > 0 ? 'p' : 'm');
	//        sprintf(ts,"%03d%c.txt",abs(mu), mu > 0? 'p' : 'm');
        // for(n=0; n < nVarPar; n++) dat >> fitpar[mu][n];
      }
      else {
        // dat.open((fpath+"0.txt").c_str());
        strcpy(ts,"0.txt");
      }
      #ifdef DEBUG_OFFSET
        cout << "--> FillFromFiles, mu = " << mu << endl;
        cout << "--> FillFromFiles, ts = " << ts << endl;
        // for(int j=0; j < nVarPar; j++) cout << j <<": "<< p[j] << endl;
      #endif
      dat.open((fpath+ts).c_str());
      if(!dat) return 1;
      dat >> v;
      if(!dat.good()) return 1;
      if(v.dim() != nVarPar) return 2;
      NewParams(v,mu);
      dat.close();
    }
    
    return 0;
  }
  
  // ============================================================
  double Param(int j) {return fitpar0[j];}
  
  // // ============================================================
  // double Error(int j) {return sqrt(Cov[j][j]);}

  #if 0
  // ============================================================
  void GetParams(double* p) {
    memcpy(p, fitpar0, nVarPar*sizeof(double));
  }
  
  // ============================================================
  void GetErrors(double* p) {
    for(int j=0; j < nVarPar; j++) p[j] = sqrt(Cov[j][j]);
  }
  #endif
  
};

#endif
