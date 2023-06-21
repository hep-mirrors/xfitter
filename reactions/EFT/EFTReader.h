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


class EFTReader {

public:
  vector<string> filename_list;
  string input_type;

  int debug = -1;
  bool abs_output = false;
  bool no_central = false;
  int rows_before_transpose = -1; // do not tranpose xsec if <= 0

  vector<string> name_EFT_param; 
  size_t num_param; // number of EFT paramters; should be less than 99
  map<string, size_t> find_EFT_param_id1; // , find_EFT_param_id0;

  size_t num_bin; // number of bins
  bool normQ = false; // todo; what do we need?

  bool scaleQ = false;
  valarray<double> final_scaling; // use valarray to simplify multiplication

  std::map<size_t, std::vector<double>* > coeff; // linear and quadratic coefficients of all EFT parameters; for fixed input
  std::map<size_t, Vec* > basis; // for mixed input

  RawVec* prvec_C = nullptr;
  vector<RawVec* > raw_basis;

  double val_EFT_param[99]; // MAX_NUM_PARAM

  /////////////////////////////////////////////////////////////////////////////
  // initialization
  EFTReader(vector<string> fname_list, string input_type_in) {
    for (string fname : fname_list)
      filename_list.push_back(fname);
    input_type = input_type_in;
  }
    
  EFTReader(vector<string> fname_list, string input_type_in, ReactionTheory* reaction) {
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
  const size_t MAX_NUM_PARAM = 99;
  void solvelq(Vec*, Vec*, RawVec*, RawVec*);
  void solvelQ(Vec*, Vec*, RawVec*, RawVec*);
  void solveNol(Vec*, Vec*, RawVec*, RawVec*);
  void lqQCoeff(vector<double>& c, int type, double val);

  void transpose(valarray<double>& xsec);
  void scaleXSec(valarray<double>& xsec);
  void calcXSecMixed(valarray<double>& xsec);
  void calcXSecFixed(valarray<double>& xsec);
  void readFixedInput();
  void readMixedInput();
};

#endif
