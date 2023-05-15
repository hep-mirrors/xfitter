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
			//--------------------------------------------

using namespace std;

class EFTReader {

public:

  // initialization
  EFTReader(vector<string> fname_list) {
    for (string fname : fname_list)
      filename_list.push_back(fname);
  }
    
  EFTReader(vector<string> fname_list, ReactionTheory* reaction) {
    for (string fname : fname_list)
      filename_list.push_back(fname);
    _reactionTheory=reaction;
  }

  // name of file  
  // string file_EFT;
  vector<string> filename_list;

  int MAX_NUM_PARAM = 100; // todo set as a global const

  int num_bin; // number of bins
  int num_param; // number of EFT paramters; should be less than 100
  std::map<int, std::vector<double>* > coeff; //a7: linear and quadratic coefficients of all val_EFT_params.
  double val_EFT_param[100]; // MAX_NUM_PARAM

  // initialization
  void setinit(vector<string> name_EFT_param);

  void setValEFT(vector<double> coe){
    // executed for each computation
    if (num_param == coe.size()) {
      for (int i=0; i<num_param; i++)
	val_EFT_param[i] = coe[i];
    } else {
      hf_errlog(23040301, "E: number of EFT parameters does not match");
    }
  };

  // calculating xsecs 
  vector<double> calcxsec();

private:

  // commons
  ReactionTheory* _reactionTheory;

};

#endif
