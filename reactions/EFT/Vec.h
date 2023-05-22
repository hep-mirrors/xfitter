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
  ///////////////////////////////////////////////////////
  int type = 0;
  vector<ingredient> ingredients;

  ///////////////////////////////////////////////////////
  Vec(int type_in);
  void addIng(RawVec* prvec, double coeff);
  void book(double val);
};

/////////////////////////////////////////////////////////////////////////////

class RawVec {
 public:
  RawVec(YAML::node node, string key);

  ///////////////////////////////////////////////////////
  int type = 4;
  string format;
  string param_name1;
  string param_name2;
  double param_val1;
  double param_val2;
  vector<string> grid_file_list;
  vector<double> ratio_list;
  vector<double> value_list; // cross sections in each bin
  double coeff = 0.0;
  ///////////////////////////////////////////////////////
  void FR2FA(vector<double> val_list_C);

  void convolute();

  void setCoeff(double val) { coeff = val; }

  void increaseCoeff(double val) { coeff += val; }

  void increaseXSecInPlace(vector<double> xsec);

  /////////////////////////////////
  // private:

};

#endif
