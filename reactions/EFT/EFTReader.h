#ifndef EFT_H
#define EFT_H

//--------------------------------------------------------------
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
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
  string inputType;

  int debug = -1;
  bool abs_output = false;
  bool no_central = false;

  vector<string> name_EFT_param; 
  int num_param; // number of EFT paramters; should be less than 99
  map<string, int> find_EFT_param_id;

  int num_bin; // number of bins

  std::map<int, std::vector<double>* > coeff; // linear and quadratic coefficients of all EFT parameters; for fixed input
  std::map<int, Vec* > basis; // for mixed input

  RawVec* prvec_C;
  vector<RawVec* > raw_basis;

  double val_EFT_param[99]; // MAX_NUM_PARAM

  /////////////////////////////////////////////////////////////////////////////
  // initialization
  EFTReader(vector<string> fname_list, string inputType_in) {
    for (string fname : fname_list)
      filename_list.push_back(fname);
    inputType = inputType_in;
  }
    
  EFTReader(vector<string> fname_list, string inputType_in, ReactionTheory* reaction) {
    for (string fname : fname_list)
      filename_list.push_back(fname);

    inputType = inputType_in;
    _reactionTheory=reaction;
  }

  void initParamName(vector<string> name_EFT_param_in);
  void read_input();
  void read_fixed_input();
  void read_mixed_input();
  void initlq();
  void initm();
  void initrvec();
  void initIter(vector<double> list_val);
  void updatervec();  
  void book();
  void setValEFT(vector<double> list_val);

  // calculating cross sections
  vector<double> calcXSecMixed();
  vector<double> calcXSecFixed();
  vector<double> calcXSec();

private:
  // commons
  ReactionTheory* _reactionTheory;
  const int MAX_NUM_PARAM = 99;
};

#endif
