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
// for fast grids
#ifdef WITH_PINEAPPL
#include "pineappl_capi.h"
#endif
#include "appl_grid/appl_grid.h"
//--------------------------------------------------------------
using namespace std;

class RawVec {
 private:
  double coeff = 0.0;
  vector<string> grid_file_list;
  vector<double> ratio_list;
  int pdg_id = 2212;
  double xi_ren = 1.0; // renom. scale
  double xi_fac = 1.0; // fact. scale
  size_t num_bin;
  bool save_grid_in_memory = true;

 public:
  int type = 4;
  string format;
  string entry;
  string param_name1;
  string param_name2;
  double param_val1;
  double param_val2;
  vector<double> value_list; // cross sections in each bin
#ifdef WITH_PINEAPPL
  vector<pineappl_grid* > pgrid_list;
#endif
  // vector<unique_ptr<appl::grid> > p_APPLgrid_list;
  vector<appl::grid* > p_APPLgrid_list;

  ///////////////////////////////////////////////////////
  RawVec (YAML::Node node, string key, size_t num_bin, string grid_dir, 
          double xi_ren_in, double xi_fac_in, bool save_grid_Q);

  /* RawVec (YAML::Node node, string key, size_t num_bin, string grid_dir); */

  /* RawVec (YAML::Node node, string key, size_t num_bin, string grid_dir,  */
  /*         double xi_ren_in, double xi_fac_in, bool save_grid_Q) { */
  /*   xi_ren = xi_ren_in; */
  /*   xi_fac = xi_fac_in; */
  /*   save_grid_in_memory = save_grid_Q; */
  /*   RawVec(node, key, num_bin, grid_dir); */
  /* } */

  void FR2FA(vector<double> val_list_C);

  void convolute();

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
  // vector<ingredient> ingredients;
  vector<ingredient*> ingredients;

  ///////////////////////////////////////////////////////
  Vec(int type_in);
  void addIng(RawVec* prvec, double coeff);
  void book(double val);
};

/////////////////////////////////////////////////////////////////////////////

#endif
