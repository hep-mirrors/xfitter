   /*
     @file EFTTerm.h
     @date 2023-03
     @author X.M. Shen <xmshen137@gmail.com>
   */
#ifndef EFT_H
#define EFT_H

//--------------------------------------------------------------
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <valarray>
#include <cstdlib>
#include <map>
#include "Vec.h"
#include "hf_errlog.h"
#include "ReactionTheory.h"
// #include "TermData.h"

using namespace std;


class EFTTerm {

public:
  string input_type;
  vector<string> filename_list;

  //
  vector<string> name_EFT_param;
  size_t num_param; // number of EFT parameters; should be less than 99
  map<string, size_t> find_EFT_param_id1; // , find_EFT_param_id0;

  //
  bool abs_output = false;
  bool no_central = false;
  double xi_ren = 1.0, xi_fac = 1.0;
  int debug = -1;

  //
  size_t num_bin; // number of bins

  //
  int rows_before_transpose = -1; // do not transpose xsec if <= 0
  bool normQ = false;
  valarray<double> binning_for_norm; // use valarray to simplify multiplication

  bool scaleQ1 = false;
  bool scaleQ2 = false;
  valarray<double> scaling1; // use valarray to simplify multiplication
  valarray<double> scaling2; // use valarray to simplify multiplication

  //
  std::map<size_t, std::vector<double>* > coeff; // linear and quadratic coefficients of all EFT parameters; for fixed input
  std::map<size_t, Vec* > basis; // for mixed input
  vector<Vec* > pvec_list;

  RawVec* prvec_C = nullptr;
  vector<RawVec* > raw_basis;

  double val_EFT_param[99]; // MAX_NUM_PARAM

  /////////////////////////////////////////////////////////////////////////////
  // initialization
  EFTTerm(vector<string> fname_list, string input_type_in) {
    for (string fname : fname_list)
      filename_list.push_back(fname);

    input_type = input_type_in;
  }

  EFTTerm(vector<string> fname_list, string input_type_in, ReactionTheory* reaction) {
    for (string fname : fname_list)
      filename_list.push_back(fname);

    input_type = input_type_in;
    _reactionTheory=reaction;
  }

  void initParamName(vector<string> name_EFT_param_in);
  void readInput();
  void initlq(size_t);
  int  initm(size_t, size_t);
  void initrvec();
  void initIter(valarray<double>& list_val);
  void updatervec();
  void book();
  void setValEFT(valarray<double>& list_val);

  // calculating cross sections
  void calcXSec(valarray<double>& xsec);

private:
  // commons
  ReactionTheory* _reactionTheory;
  const size_t MAX_NUM_PARAM = 99; // todo
  void solvelq(Vec*, Vec*, RawVec*, RawVec*);
  void solvelQ(Vec*, Vec*, RawVec*, RawVec*);
  void solveNol(Vec*, Vec*, RawVec*, RawVec*);
  void lqQCoeff(vector<double>& c, int type, double val);

  void transpose(valarray<double>& xsec);
  void scaleXSec1(valarray<double>& xsec);
  void scaleXSec2(valarray<double>& xsec);
  void normXSec(valarray<double>& xsec);
  void calcXSecMixed(valarray<double>& xsec);
  void calcXSecFixed(valarray<double>& xsec);
  void readFixedInput();
  void readMixedInput();
};

#endif
