// #include <vector>
// #include <string>
// #include <cstring>
// #include <iostream>
// #include <fstream>
#include "yaml-cpp/yaml.h"
#include "Vec.h"

//--------------------------------------------------------------
using namespace std;

Vec::Vec(int type_in) {
  type = type_in;
  if (type > 3 || type < 1)
    std::cout << "valid type = 1,2,3 for l,q,m" << std::endl;
  // hf_errlog(23051602, "valid type = 1,2,3 for l,q,m");
}

void Vec::addIng(RawVec* prvec, double coeff) {
  // ingredient* ping = new ingredient { prvec, coeff };
  bool findQ = false;

  if (type == 3) {
    for (auto ing: ingredients) {
      if (ing->prvec == prvec) {
	ing->coeff += coeff;
	findQ = true;
	break;
      }
    }
  }
  if (findQ == false){
    ingredient* p = new ingredient();
    p->prvec = prvec;
    p->coeff = coeff;
    ingredients.push_back(p);
  }
}
///////////////////////////////////////////////////////
void Vec::book(double val) {
  for (auto ing: ingredients)
    ing->prvec->increaseCoeff(ing->coeff * val);
}


/////////////////////////////////////////////////////////////////////////////
// RawVec
/////////////////////////////////////////////////////////////////////////////
RawVec::RawVec (YAML::Node node, string key) {
  // RawVec::RawVec (YAML::const_iterator node, string key) {
  // key: tag for the current entry; only used for issuing errors

  // read type
  if (node["type"]) {
    string typeS =  node["type"].as<string>();
    if      (typeS == "C")
      type = 0;
    else if (typeS == "l")
      type = 1;
    else if (typeS == "L")
      type = -1;
    else if (typeS == "q")
      type = 2;
    else if (typeS == "Q")
      type = -2;
    else if (typeS == "m")
      type = 3;
    else if (typeS == "M")
      type = -3;
    else
      cout << "invalid type" << endl;
  }
  else {
    cout << "Error: entry type not given: " << key << endl;
  }
  ///////////////////////////////////////////////////////
  // read input
  if (node["format"]) 
    format = node["format"].as<string>();
  else
    cout << "Error: entry format not given: " << key << endl;

  if (node["xsec"]) {
    if (format == "FR") {
      ratio_list = node["xsec"].as<vector<double> >();
      for (size_t i=0; i<ratio_list.size(); i++)
	value_list.push_back(0.0);
    }
    else if (format == "FA")
      value_list = node["xsec"].as<vector<double> >();
    else if (format == "PineAPPL" || format == "fastNLO" || format == "APPLgrid")
      grid_file_list = node["xsec"].as<vector<string> >();
    else
      cout << "Error: grid format not support" << endl;
  }
  else
    cout << "Error: values(grids) for fixed input(mixed) not given: " << key  << endl;

  ///////////////////////////////////////////////////////      
  // read EFT parameters
  if (type == -3 || type == 3) {
    // names of parameters
    if (node["param"]) {
      vector<string> params = node["param"].as<vector<string> >();
      if (params.size() != 2) 
	cout << "Error: number of parameters" << endl;
      else {
	param_name1 = params[0];
	param_name2 = params[1];
      }
    } 
    else
      cout << "Error: param not found" << endl;

    if (type == -3) {
      // value of parameters
      if (node["param_value"]) {
	vector<double> param_vals = node["param_value"].as<vector<double> >();

	if (param_vals.size() != 2) 
	  cout << "Error: number of parameters" << endl;
	else {
	  param_val1 = param_vals[0];
	  param_val2 = param_vals[1];
	}
      } 
      else
	cout << "Error: param not found" << endl;
    }
  } // end of 3,-3
  else if (type != 0) {
    // names of parameter
    if (node["param"])
      param_name1 = node["param"].as<string>();
    else
      cout << "Error: param not found" << endl;

    if (type < 0) {
      // value of parameter
      if (node["param_value"])
	param_val1 = node["param_value"].as<double>();
      else
	cout << "Error: param not found" << endl;
    }
  }
} // end of constructor


/////////////////////////////////

void RawVec::FR2FA(vector<double> val_list_C) {
  // ratio -> absolute value
  assert(format == "FR");

  if (value_list.size() != val_list_C.size()) 
    cout << "Error: size does not match:" << value_list.size() << ", " << val_list_C.size()  << endl;
  else {
    for (size_t i=0; i<val_list_C.size(); i++)
      value_list[i] = ratio_list[i] * val_list_C[i];
  }
}

void RawVec::convolute() {
  cout << "Error: convoute() not implemented" << endl;
} // todo

/////////////////////////////////
void RawVec::increaseXSecInPlace(valarray<double> xsec) {
  if (xsec.size() != value_list.size())
    cout << "Error: size does not match:" << xsec.size() << ", " << value_list.size() << endl;
  else {
    for (size_t i=0; i<value_list.size(); i++)
      xsec[i] += value_list[i] * coeff;
  }
}
/////////////////////////////////
