#ifndef EFT_H
#define EFT_H

//--------------------------------------------------------------
// C++ wraper for using EFT interface 

#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "hf_errlog.h"
#include "ReactionTheory.h"
#include "TermData.h"
#include <map>
#include "Vec.h"
			//--------------------------------------------

using namespace std;


class EFTReader {

public:

  // initialization
  EFTReader(vector<string> fname_list, int debug_in) {
    for (string fname : fname_list)
      filename_list.push_back(fname);
    debug = debug_in;
  }
    
  EFTReader(vector<string> fname_list, int debug_in, ReactionTheory* reaction) {
    for (string fname : fname_list)
      filename_list.push_back(fname);

    debug = debug_in;
    _reactionTheory=reaction;
  }

  vector<string> filename_list;
  int debug = -1;

  int MAX_NUM_PARAM = 100; // todo set as a global const

  int num_bin; // number of bins
  int num_param; // number of EFT paramters; should be less than 100
  vector<string> name_EFT_param; 
  double val_EFT_param[100]; // MAX_NUM_PARAM
  std::map<int, std::vector<double>* > coeff; // linear and quadratic coefficients of all EFT parameters

  // initialization
  void init(vector<string> name_EFT_param);

  void setValEFT(vector<double> coe) {
    // executed for each computation
    if (num_param == coe.size()) {
      for (int i=0; i<num_param; i++)
	val_EFT_param[i] = coe[i];
    } else {
      hf_errlog(23040301, "E: number of EFT parameters does not match");
    }

    if (debug > 2) {
      std::cout << "=======================================================" << std::endl;
      std::cout << "EFTReader.setValEFT" << std::endl;
      for (int i=0; i<num_param; i++) {
	std::cout << name_EFT_param[i] << "=" <<  val_EFT_param[i] << std::endl;
      }
    }
  };

  // calculating xsecs 
  vector<double> calcXSec();

private:
  // commons
  ReactionTheory* _reactionTheory;

};

#endif
