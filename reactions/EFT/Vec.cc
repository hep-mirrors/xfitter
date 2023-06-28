// #include <vector>
// #include <string>
// #include <cstring>
// #include <iostream>
// #include <fstream>
#include "yaml-cpp/yaml.h"
#include "Vec.h"

// from ReactionAPPLgrid
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include "xfitter_cpp_base.h"
#include <memory>
#include "BaseEvolution.h"

#include "appl_grid/appl_grid.h"
#include "TermData.h" // needed by the PDF wrappers
// from Toni's script:
#include "pineappl_capi.h"
#include <sstream>


//--------------------------------------------------------------
using namespace std;

Vec::Vec(int type_in) {
  type = type_in;
  if (type > 3 || type < 1) {
    hf_errlog(23051602, "S: valid type = 1,2,3 for l,q,m");
  }
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
  for (auto ing: ingredients) {
    ing->prvec->increaseCoeff(ing->coeff * val);
  }
}


/////////////////////////////////////////////////////////////////////////////
// RawVec
/////////////////////////////////////////////////////////////////////////////
//RawVec::RawVec (YAML::Node node, string key, size_t num_bin_in, string grid_dir) {
RawVec::RawVec (YAML::Node node, string key, size_t num_bin_in, string grid_dir, 
          double xi_ren_in, double xi_fac_in, bool save_grid_Q) {
  // key: tag for the current entry; only used for issuing errors

  xi_ren = xi_ren_in;
  xi_fac = xi_fac_in;
  save_grid_in_memory = save_grid_Q;
  entry = key;
  num_bin = num_bin_in;
  ///////////////////////////////////////////////////////
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
      hf_errlog(23061505, "S: invalid type for entry " + key);
  }
  else {
    hf_errlog(23061505, "S: type not given for entry " + key);
  }
  ///////////////////////////////////////////////////////
  // read xiF, xiR, pdg_id
  if (node["xi_ren"]) 
    setScaleRen(node["xi_ren"].as<double>());
  if (node["xi_fac"]) 
    setScaleFac(node["xi_fac"].as<double>());
  if (node["pdg_id"]) 
    setPDGId(node["pdg_id"].as<int>());


  ///////////////////////////////////////////////////////
  // read xsec
  if (node["format"]) {
    format = node["format"].as<string>();
  }
  else
    hf_errlog(23061504, "S: format not given for entry " + key);

  if (node["xsec"]) {
    if (format == "FR") {
      ratio_list = node["xsec"].as<vector<double> >();

      if (ratio_list.size() != num_bin) {
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << ratio_list.size() << " v.s. " <<  num_bin  << endl;
	hf_errlog(23061501, "S: length of xsec does not match info:num_bin for entry " + key);
      }

      for (size_t i=0; i<num_bin; i++)
	value_list.push_back(0.0);
    }
    else if (format == "FA") {
      value_list = node["xsec"].as<vector<double> >();

      if (value_list.size() != num_bin) {
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << value_list.size() << " v.s. " <<  num_bin  << endl;
	hf_errlog(23061521, "S: length of xsec does not match info:num_bin for entry " + key);
      }
    }
    else if (format == "PineAPPL" || format == "fastNLO" || format == "APPLgrid") {
      for (size_t i=0; i<num_bin; i++)
	value_list.push_back(0.0);

      for (string filepath: node["xsec"].as<vector<string> >())
	grid_file_list.push_back(grid_dir + "/" + filepath);
    }
    else
      hf_errlog(23061503, "S: grid format not support for entry " + key);
  }
  else
    hf_errlog(23061502, "S: cross section(grids) for fixed input(mixed) not given for entry " + key);

  ///////////////////////////////////////////////////////      
  // read grids
  // todo: check if num_bin matches
  if ( save_grid_in_memory ) {
    // todo check format type
    for (string grid_file_name: grid_file_list) {
      if (format == "PineAPPL") {
	pineappl_grid* g = pineappl_grid_read(grid_file_name.c_str());
	pgrid_list.push_back(g);
	hf_errlog(23061202, "I: read PineAPPL grid from " + grid_file_name);
      }
      else if (format == "APPLgrid") {
	
      }
    }
  }  
  ///////////////////////////////////////////////////////      
  // read EFT parameters
  if (type == -3 || type == 3) {
    // names of parameters
    if (node["param"]) {
      vector<string> params = node["param"].as<vector<string> >();
      if (params.size() != 2) 
	hf_errlog(23061510, "S: check `param` for entry:" +entry);
      else {
	param_name1 = params[0];
	param_name2 = params[1];
      }
    } 
    else
      hf_errlog(23061511, "S: `param` not found for entry:" +entry);

    if (type == -3) {
      // value of parameters
      if (node["param_value"]) {
	vector<double> param_vals = node["param_value"].as<vector<double> >();

	if (param_vals.size() != 2) 
	  hf_errlog(23061513, "S: `param_value` for entry:" +entry);
	else {
	  param_val1 = param_vals[0];
	  param_val2 = param_vals[1];
	}
      } 
      else
	hf_errlog(23061511, "S: `param_value` not found for entry:" +entry);
    }
  } // end of 3,-3
  else if (type != 0) {
    // names of parameter
    if (node["param"])
      param_name1 = node["param"].as<string>();
    else
      hf_errlog(23061511, "S: `param` not found for entry:" +entry);

    if (type < 0) {
      // value of parameter
      if (node["param_value"])
	param_val1 = node["param_value"].as<double>();
      else
	hf_errlog(23061511, "S: `param_value` not found for entry:" +entry);
    }
  }


} // end of constructor


/////////////////////////////////

void RawVec::FR2FA(vector<double> val_list_C) {
  // ratio -> absolute value
  assert(format == "FR");

  if (value_list.size() != val_list_C.size()) 
    hf_errlog(23061515, "S: size does not match");
  else {
    for (size_t i=0; i<val_list_C.size(); i++)
      value_list[i] = ratio_list[i] * val_list_C[i];
  }
}

void RawVec::convolute() {
  if (format != "PineAPPL") {
    hf_errlog(23061201, "S: Grids other than PineAPPL are not supported yet");
  }
  /////////////////////////////////////////////////////////////////////////////
  // for PineAPPL
  /////////////////////////////////////////////////////////////////////////////
  if (format == "PineAPPL") {
    // read the grids
    // todo: follow Toni's code to deal with exceptions
    if ( ! save_grid_in_memory ) {
      for (string grid_file_name: grid_file_list) {
	pineappl_grid* g = pineappl_grid_read(grid_file_name.c_str());
	pgrid_list.push_back(g);
	hf_errlog(23061202, "I: read PineAPPL grid from " + grid_file_name);
      }
    }
    ////////////////////////////////////////////
    // convolute
    ////////////////////////////////////////////
    // PDFs can be called with wrappers within xFitter
    // Here we follow Toni's strategy to deal with PDFs
    auto xfx = [](int32_t id_in, double x, double q2, void *state) {
      double pdfs[13];
      int32_t id = id_in==21 ? 6 : id_in+6;
      pdf_xfxq_wrapper_(x, sqrt(q2), pdfs);
      return pdfs[id];
    };
    auto alphas = [](double q2, void *state) {
      return alphas_wrapper_(sqrt(q2));
    };

    //See function specification in deps/pineappl/include/pineappl_capi/pineappl_capi.h
    int shift_bins = 0;

    for (auto pgrid: pgrid_list) {
      pineappl_grid_convolute_with_one(
                         pgrid, pdg_id,
                         xfx, alphas,
                         nullptr, //"state" provided to wrappers, redundant in xFitter
                         nullptr, // order mask
                         nullptr, // lumi mask
                         xi_ren, xi_fac,
                         value_list.data() + shift_bins);

      shift_bins += pineappl_grid_bin_count(pgrid);
    } // end of loop over all grid files

    // free the grids
    if (! save_grid_in_memory) {
      for (auto p: pgrid_list)
	pineappl_grid_delete(p);
      pgrid_list.clear();
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  // for APPLgrid
  if (format == "APPLgrid") {
  }
}

/////////////////////////////////
void RawVec::increaseXSecInPlace(valarray<double>& xsec) {
  if (xsec.size() != value_list.size())
    hf_errlog(23061516, "S: size does not match");
  else {
    for (size_t i=0; i<value_list.size(); i++) {
      xsec[i] += value_list[i] * coeff;
    }
  }
}
/////////////////////////////////
