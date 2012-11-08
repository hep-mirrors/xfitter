/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
  \date 2011--2012
  \copyright Creative Commons license CC-BY-NC 3.0
_____________________________________________________________*/

#ifndef _DECOR_H_
#define _DECOR_H_

#include "matrix.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

// #define DBGT_ON
#include "dbgtools.h"
// #define DEBUG_OFFSET

//ooooooooooooooooooooooooooooooooooooooooo
class Decor_t {

  int N;
  SqMatrix_t Cov;
  SqMatrix_t EigVecs; // = matrix diagonalizing Cov;
  Vector EigVals;
  Vector EigShifts;
  // ==>  EigVecs[i][m] * Cov[i][j] * EigVecs[j][n] = delta(m,n)*EigVals[n]
  // --- shift of a[i] by one std. dev. in the direction j-th decorrelated error = EigVecs[i][j]*sqrt(EigVals[j])
  bool AutoDiag;
  double DeltaChiSqr;
  
public:
  //===========================
  Decor_t() : N(0), AutoDiag(false), DeltaChiSqr(1.0) {}
  
  //===========================
  void Diag() {
    DBG_SHOW(Cov)
    EigVecs = Cov;
    EigVals = Vector(0, N);
    EigVecs.Diagonalize(EigVals);
    // cout << EigVals << endl;
    EigShifts = Vector(0, N);
    for(int j=0; j < N; j++) EigShifts[j] = sqrt(DeltaChiSqr*EigVals[j]);
    // SetDeltaChiSq(DeltaChiSqr);
  }

  //===========================
  void SetCov(const SqMatrix_t& C) {
    Cov = C;
    N = C.len;
    if(AutoDiag) Diag();
  }

  //===========================
  void SetCov(int n, double *v) {
    Cov = SqMatrix_t(0, N=n);
    Cov.Set(v);
    if(AutoDiag) Diag();
  }

  //===========================
  void SetDeltaChiSq(double dc2) {
    DeltaChiSqr = dc2;
  }

  //===========================
  int GetN() const {return N;}

  //===========================
  int SetCov(const char* fn);

  //===========================
  void Sort(bool asc=1);

  //===========================
  double GetEigShift(int iPropErr) const {
    return EigShifts[iPropErr];
  }

  //===========================
  double VarError(int iVar, int iPropErr) const {
    // return EigVecs[iVar][iPropErr]*sqrt(DeltaChiSqr*EigVals[iPropErr]);
    return EigVecs[iVar][iPropErr]*EigShifts[iPropErr];
  }

};

#endif
