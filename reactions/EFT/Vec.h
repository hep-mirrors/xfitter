   /*
     @file Vec.h
     @date 2023-03
     @author X.M. Shen <xmshen137@gmail.com>
   */
#ifndef VECH
#define VECH

//--------------------------------------------------------------
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <valarray>
#include <cassert>
#include "yaml-cpp/yaml.h"
#include <utility> // For std::pair
// for fast grids
#ifdef WITH_PINEAPPL
#include "pineappl_capi.h"
#endif

#include "appl_grid/appl_grid.h"

// #include "fastnlotk/fastNLOReader.h" // need the fastnlotk directory under the current dir.
//--------------------------------------------------------------
using namespace std;

// format
const string FRstr = "ratio"; // also needed by EFTTerm.cc
const string FAstr = "xsection";
const string PineAPPL = "PineAPPL";
const string APPLgrid = "APPLgrid";

// if >= 0; do not need param_value
const int typeC = 0;
const int typel = 1;
const int typeq = 2;
const int typem = 11; // arbitrary
const int typeMonoMin = 3;
const int typeMonomial = 137; // arbitrary
// if < 0; need param_value
const int typeL = -1;
const int typeQ = -2;
const int typeM = -11;

class RawVec {
 private:
  double coeff = 0.0;
  vector<string> grid_file_list;
  vector<double> ratio_list;
  int pdg_id = 2212;
  double xi_ren = 1.0; // renom. scale
  double xi_fac = 1.0; // fact. scale
  size_t num_bin;
  // size_t* p_num_bin = &num_bin;
  bool save_grid_in_memory = true;
  void initReadEFTParam(YAML::Node node);
 public:
  int type = 4;
  string format;
  string entry;
  string param_name1;
  string param_name2;
  double param_val1;
  double param_val2;
  vector<double> value_list; // cross sections in each bin
  vector<string> param_name_list; // for monomials
  vector<int> power_list; // for monomials
  // PineAPPL grids
#ifdef WITH_PINEAPPL
  vector<pineappl_grid* > pgrid_list;
#endif
  // APPLgrid grids
  // vector<unique_ptr<appl::grid> > p_APPLgrid_list;
  vector<appl::grid* > p_APPLgrid_list;
  // fastNLO grids
  // todo
  // vector< *> p_fastNLO_list;
  ///////////////////////////////////////////////////////
  RawVec(YAML::Node node, string key, size_t & num_bin, string grid_dir, 
         double xi_ren_in, double xi_fac_in, bool save_grid_Q);

  /* RawVec (YAML::Node node, string key, size_t num_bin, string grid_dir); */

  void FR2FA(vector<double> val_list_C);

  void convolute();
  void convolute_PineAPPL();
  void convolute_APPLgrid();
  // void convolute_fastNLO();

  void setCoeff(double val) { coeff = val; }

  void increaseCoeff(double val) { coeff += val; }

  void increaseXSecInPlace(valarray<double>& xsec);

  void setScaleRen(double xi_ren_in) {xi_ren = xi_ren_in; }

  void setScaleFac(double xi_fac_in) {xi_fac = xi_fac_in; }

  void setPDGId(int id) {pdg_id = id; }
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
  vector<std::pair<int, int> > param_id_power_list;
  vector<ingredient*> ingredients;

  ///////////////////////////////////////////////////////
  Vec(int type_in);
  void addIng(RawVec* prvec, double coeff);
  void book(double val);
  // void bookMonomial(valarray<double>& val);
  void bookMonomial(double val[99]);
  void addParamPower(int param_id, int power);
};


/////////////////////////////////////////////////////////////////////////////

#endif
