#ifndef VECH
#define VECH

//--------------------------------------------------------------
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include <cassert>

//--------------------------------------------------------------
using namespace std;

class RawVec {
 private:
  int type = 4;
  string format;
  double coeff = 0.0;
  vector<string> grid_file_list;
  vector<double> ratio_list;

 public:
  string param_name1;
  string param_name2;
  double param_val1;
  double param_val2;
  vector<double> value_list; // cross sections in each bin

  ///////////////////////////////////////////////////////
  RawVec (YAML::Node node, string key);

  void FR2FA(vector<double> val_list_C);

  void convolute();

  void setCoeff(double val) { coeff = val; }

  void increaseCoeff(double val) { coeff += val; }

  void increaseXSecInPlace(vector<double> xsec);

};

/////////////////////////////////////////////////////////////////////////////

struct ingredient{
  RawVec* prvec;
  double coeff=0.0;
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

#endif
