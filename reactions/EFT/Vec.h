#ifndef VECH
#define VECH

//--------------------------------------------------------------
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>

//--------------------------------------------------------------
using namespace std;

struct ingredient{
  RawVec* prvec;
  double coeff=0;
};

/////////////////////////////////////////////////////////////////////////////
class Vec {
 public:
  // constructor
  Vec(int type_in) {
    type = type_in;
    if (type > 3 || type < 1)
      std::cout << "valid type = 1,2,3 for l,q,m" << std::endl;
      // hf_errlog(23051602, "valid type = 1,2,3 for l,q,m");
  }
  ///////////////////////////////////////////////////////
  int type = 0;
  vector<ingredient> ingredients;
  ///////////////////////////////////////////////////////
  void addIng(RawVec* prvec, double coeff) {
    // ingredient* ping = new ingredient { prvec, coeff };
    bool findQ = false;

    if (self.type == 3) {
      for (auto ing: ingredients) {
	if (ing.prvec == prvec) {
	  ing.coeff += coeff;
	  findQ = true;
	  break;
	}
      }
    }
    if (findQ == false)
      ingredients.push_back(ingredient{ prvec, coeff });
  }
  ///////////////////////////////////////////////////////
  void book(double val) {
    for (auto ing: ingredients)
      ing.prvec->increase_coeff(ing.coeff * val);
  }
  ///////////////////////////////////////////////////////
  //private:
};

/////////////////////////////////////////////////////////////////////////////

class RawVec {
 public:
  RawVec() {
  }
  ///////////////////////////////////////////////////////
  int type = 4;
  string format;
  double param_val1;
  double param_val2;
  string param_name1;
  string param_name2;
  vector<string> grid_file_list;
  vector<double> value_list;
  double coeff = 0.0;
  ///////////////////////////////////////////////////////
  void increaseXSecInPlace() {
  }

  void setCoeff(double val) { coeff = val; }
  void increaseCoeff(double val) { coeff += val; }
  void convolute() {}
  void FR2FA() {}

  // private:

};

#endif
