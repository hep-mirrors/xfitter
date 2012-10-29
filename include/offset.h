/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _OFFSET_H_
#define _OFFSET_H_

// #define DEBUG_OFFSET

#include <cstdio>
#include <fstream>
#include "matrix.h"

// ooooooooooooooooooooooooooo
class OffsetCalc_t {
public:
  // typedef vector <double> rvec_t;
  typedef double*  RealPtr_t;
  RealPtr_t *fitpar, *fitpar_0;
  // vector <rvec_t> DPC;
  int nCosy;
  int nVarPar;
  SqMatrix_t Cov, CosyCov, StatCov;
  
  // ============================================================
  OffsetCalc_t(int n_Cosy, int n_VarPar) : nCosy(n_Cosy), nVarPar(n_VarPar) {
    // --- create matrix for shifted parameters
    int nC = 2*nCosy+1;
    fitpar_0 = new RealPtr_t[nC];
    for(int iq = 0; iq < nC; iq++) {
      fitpar_0[iq] = new double[nVarPar];
      // memset(fitpar_0[iq],0,nCosy*sizeof(double));
    }
    fitpar = fitpar_0 + nCosy;
    
    // --- create covariance matrices
    Cov = SqMatrix_t(0, nVarPar);
    CosyCov = SqMatrix_t(0, nVarPar);
    StatCov = SqMatrix_t(0, nVarPar);
    #ifdef DEBUG_OFFSET
      cout << "--> OffsetCalc_t: " << "nCosy = "<< nCosy <<", nVarPar = "<< nVarPar << endl;
    #endif
  }
  
  // ============================================================
  ~OffsetCalc_t() {
    int nC = 2*nCosy+1;
    for(int iq = 0; iq < nC; iq++) delete[] fitpar_0[iq];
    delete[] fitpar_0;
  }
  
  // ============================================================
  void NewParams(double* p, int mu) {
    memcpy(fitpar[mu], p, nVarPar*sizeof(double));
    #ifdef DEBUG_OFFSET
      cout << "--> NewParams, mu = " << mu << endl;
      for(int j=0; j < nVarPar; j++) cout << j <<": "<< fitpar[mu][j] << endl;
    #endif
  }
  
  // ============================================================
  void SetCovStat(double* p) {
    StatCov.Set(p);
    // StatCov.Transpose();
    #ifdef DEBUG_OFFSET
      cout << "--> SetCovStat" << endl;
      cout << StatCov << endl;
    #endif
  }
  
  // // ============================================================
  // void FillDPC(double* p) {
    // for(int j=0; j < nVarPar; j++)
      // for(int mu=1; mu <= nCosy; mu++) DPC[j][mu] = (fitpar[mu][j] - fitpar[-mu][j])/2;
  // }
  
  // ============================================================
  double DPC(int j, int mu) {
    return (fitpar[mu][j] - fitpar[-mu][j])/2;
  }
  
  // ============================================================
  void FillCov() {
    for(int j=0; j < nVarPar; j++) {
      for(int k=0; k < nVarPar; k++) {
        CosyCov[j][k] = 0;
        for(int mu=1; mu <= nCosy; mu++) CosyCov[j][k] += DPC(j,mu)*DPC(k,mu);
      }
    }
    Cov = StatCov + CosyCov;
    #ifdef DEBUG_OFFSET
      cout << "--> FillCov" << endl;
      cout << "StatCov:\n" << StatCov << endl;
      cout << "CosyCov:\n" << CosyCov << endl;
      cout << "Cov:\n" << Cov << endl;
    #endif
  }
  
  // ============================================================
  void GetCov(double* p) {
    FillCov();
    Cov.Get(p);
  }
  
  // ============================================================
  int FillFromFiles(const string& outdir) {
    string fpath(outdir);
    if(*(fpath.end()-1) != '/') fpath.append("/");
    ifstream dat;
    int n;
    
    #ifdef DEBUG_OFFSET
      cout << "FillFromFiles opening: " << (fpath+"statcov_0.txt") << endl;
    #endif
    dat.open((fpath+"statcov_0.txt").c_str());
    if(!dat) return 1;
    dat >> n;
    #ifdef DEBUG_OFFSET
      cout << "FillFromFiles n: " << n <<" = "<< nVarPar << endl;
    #endif
    if(n != nVarPar) {dat.close(); return 2;}
    double v[n=nVarPar*nVarPar], *vp = v;
    while(n--) dat >> *vp++;
    dat.close();
    SetCovStat(v);
    
    fpath.append("params_");
    char ts[16];
    for(int mu = -nCosy; mu <= nCosy; mu++) {
      if(mu) {
        sprintf(ts,"%03d%c.txt",abs(mu), mu > 0? 'p' : 'm');
        dat.open((fpath+ts).c_str());
      }
      else dat.open((fpath+"0.txt").c_str());
      if(!dat) return 1;
      dat >> n;
      #ifdef DEBUG_OFFSET
        cout << "FillFromFiles n: " << n << endl;
      #endif
      if(n != nVarPar) {dat.close(); return 2;}
      for(n=0; n < nVarPar; n++) dat >> fitpar[mu][n];
      dat.close();
    }
    
    FillCov();
    // Cov.Get(p);
    return 0;
  }
  
  // ============================================================
  void GetParams(double* p) {
    memcpy(p, fitpar[0], nVarPar*sizeof(double));
  }
  
  // ============================================================
  void GetErrors(double* p) {
    for(int j=0; j < nVarPar; j++) p[j] = sqrt(Cov[j][j]);
  }
  
};

#endif
