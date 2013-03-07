#ifndef _COV_MATRIX_H_
#define _COV_MATRIX_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

// #define DBGT_ON
#include "dbgtools.h"
#include "Xstring.h"
#include "TextReader.h"

#include "tnt_math_utils.h"
#include "tnt_matrix.h"
#include "tnt_linalg.h"

using namespace std;
using namespace TNT;

// #define VECTOR Vector<double>
// #define MATRIX Matrix<double>
typedef Vector<double> VECTOR;
typedef Matrix<double> MATRIX;

/**
  Covariance matrix
*/
// ooooooooooooooooooooooooooooooooooooooo
// template <typename T> class CovMatrix_t : public Matrix<T> {
class CovMatrix_t : public MATRIX {

  // --- DATA
  // ____________

  // int M_verb;
public:
    VECTOR EigVals;
    VECTOR EigShifts;
    MATRIX EigVecs;
  
  // --- METHODS
  // ______________
  
private:
  /**
    Read from a file. Start reading from a line matching glob pattern \c start_pat.
    The last token of this line must be the dimension of the matrix to read.
    \n
    Return error_code
  */
  // ================================================
  int read_file(const string& fname, const string& start_pat, bool LT, bool foe=true) {
    cout << "Reading Covariance Matrix from '"<< fname <<"'"<< endl;
    TextReader_t Input(fname);
    Input.SetFailOnError(foe);
    Input.Open();
    // if(Input.ErrorCode()) return Input.ErrorCode();
    if(Input.ErrorCode()) return 1;
    Xstring buf = Input.SkipUntil(start_pat.c_str());
    if(Input.ErrorCode()) return 2;

    // Input.SetFailOnError();
    StringList toks = buf.Split(" ");
    int N = toks[toks.size()-1].GetInt();
    newsize(N, N);
    double **A = *this;
    if(LT) for(int ir = 0; ir < N; ir++) {
      for(int ic = 0; ic < ir; ic++) {
        if(!(Input.Istream() >> A[ir][ic])) break;
        A[ic][ir] = A[ir][ic];
      }
      if(!(Input.Istream() >> A[ir][ir])) break;
    }
    else for(int ir = 0; ir < N; ir++) for(int ic = 0; ic < N; ic++) 
            if(!(Input.Istream() >> A[ir][ic])) break;
    return Input.Istream().good() ? 0 : 3;
  }
  

public:

  // =========================================
  CovMatrix_t& operator=(const MATRIX& A) {
    MATRIX::operator=(A);
    // EigVals   = A.EigVals;
    // EigShifts = A.EigShifts;
    // EigVecs   = A.EigVecs;
    return *this;
  }

  // =========================================
  CovMatrix_t& operator=(const CovMatrix_t& A) {
    MATRIX::operator=(A);
    EigVals   = A.EigVals;
    EigShifts = A.EigShifts;
    EigVecs   = A.EigVecs;
    return *this;
  }

  /**
    Read from a file. Start reading from a line matching glob pattern \c start_pat.
    The last token of this line must be the dimension of the matrix to read.
    \n
    Return error_code
  */
  // ================================================
  int Read(const string& fname, const string& start_pat="*") {
    // SET COVARIANCE     7
    // Covariance 7
    return read_file(fname, start_pat, false);
  }
  
  // ================================================
  int ReadL(const string& fname, const string& start_pat="*") {
    return read_file(fname, start_pat, true);
  }
  
  // ============================
  void ShowL(ostream& ostr=cout, const string& title="", bool correl_=false) {
    // cout << "--------------------------------------" << endl;
    int wid = ostr.width();
    int prec = ostr.precision();
    bool is_sci = ostr.flags() & ios::scientific;
    if(is_sci) wid = prec+8;
    else if(correl_) wid = prec+3;
    int N = num_rows();
    double **A = *this;
    if(!title.empty()) ostr << setw(0) << title << "  ";
    ostr << N << endl;
    DBG_SHOW(correl_)
    for(int ir = 0; ir < N; ir++) {
      for(int ic = 0; ic <= ir; ic++) {
        ostr << " ";
        // if(is_sci) ostr << setw(prec+8) << right;
        // else 
        ostr << right << setw(wid);
        ostr << (correl_ ? A[ir][ic]/sqrt(A[ir][ir]*A[ic][ic]) : A[ir][ic]);
      }
      ostr << endl;
    }
    ostr << endl;
    // ostr.flags(orgflags);
    // ostr.precision(orgprec);
  }
  
  // ==================================
  void deCorrelate(double DeltaChiSqr=1) {
    Linear_Algebra::Eigenvalue<double> eigen(*this);
    eigen.getRealEigenvalues(EigVals);
    eigen.getV(EigVecs);
    int N = num_rows();
    EigShifts.newsize(N);
    for(int j=0; j < N; j++) EigShifts[j] = sqrt(DeltaChiSqr*EigVals[j]);
  }
  
  //===========================
  double VarError(int iVar, int iPropErr) const {
    // return EigVecs[iVar][iPropErr]*sqrt(DeltaChiSqr*EigVals[iPropErr]);
    return EigVecs[iVar][iPropErr]*EigShifts[iPropErr];
  }
  
  // ==================================
  void FillGamma(double* fortranArray, int nrows, int ncols, int nr_) {
    DBG_SHOW(nr_)
    // DBG_SHOW(nc_)
    // --- my Gamma is transposed wrt. HFitter -- thus FTN vc C is just OK
    MATRIX Gamma; //(nrows,ncols,fortranArray);
    Gamma.Mount(ncols,nrows,fortranArray, nr_);
    // Gamma.newsize(ncols,nrows);
    // for(int i=0; i < nrows; i++) for(int j=0; j < ncols; j++) Gamma[j][i] = fortranArray[j*nr_ + i];
    // cout << scientific;
    cout << "Gamma  " << Gamma << endl;
  }
  
};

#endif
