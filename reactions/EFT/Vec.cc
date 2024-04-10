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
#ifdef WITH_PINEAPPL
#include "pineappl_capi.h"
#endif
#include <sstream>


//--------------------------------------------------------------
using namespace std;

Vec::Vec(int type_in) {
  type = type_in;
  if (type > 3 || type < 1) {
    hf_errlog(23051602, "S: EFT: valid type = 1,2,3 for l,q,m");
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
  // read xsec (numbers or filenames of grids)
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
    } // FR
    else if (format == "FA") {
      value_list = node["xsec"].as<vector<double> >();

      if (value_list.size() != num_bin) {
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << value_list.size() << " v.s. " <<  num_bin  << endl;
	hf_errlog(23061521, "S: length of xsec does not match info:num_bin for entry " + key);
      }
    } // FA
    else if (format == "PineAPPL" || format == "fastNLO" || format == "APPLgrid") {
      for (size_t i=0; i<num_bin; i++)
	value_list.push_back(0.0);

      for (string filename: node["xsec"].as<vector<string> >()) {
	std::string filepath = grid_dir + "/" + filename;
	std::ifstream file_stream(filepath);
	if (file_stream.good()) {
	  grid_file_list.push_back(filepath);
	} else {
	  std::cout << "grid file does not exist:"  << std::endl
		    << filepath << std::endl
		    << grid_dir << std::endl
		    << filename << std::endl
		    << std::endl;
	  hf_errlog(23091201, "S: grid file does not exist.");
	}
      }
    } // various grids
    else {
      hf_errlog(23061503, "S: grid format not support for entry " + key);
    }
  }
  else {
    hf_errlog(23061502, "S: cross section(grids) for fixed input(mixed) not given for entry " + key);
  }
  // end of reading cross section
  ///////////////////////////////////////////////////////      
  // read grids
  // todo: check if num_bin matches
  if ( save_grid_in_memory ) {
    for (string grid_file_name: grid_file_list) {

      if (format == "APPLgrid") {
	// XM: simply mimic ReactionAPPLgrid.cc without fully understand the code
	// not tested!
	appl::grid* g = new appl::grid(grid_file_name);
	g->trim();
	p_APPLgrid_list.push_back(g);
	hf_errlog(24040902, "I: read APPLgrid from " + grid_file_name);
      }
      else if (format == "PineAPPL") {
#ifdef WITH_PINEAPPL
	pineappl_grid* g = pineappl_grid_read(grid_file_name.c_str());
	pgrid_list.push_back(g);
	hf_errlog(23061202, "I: read PineAPPL grid from " + grid_file_name);
#else
	hf_errlog(24040901, "S: PineAPPL support is not installed");	
#endif
      }
      else if (format == "fastNLO") {
	hf_errlog(24040903, "S: EFT reaction: fastNLO support not realized");		
      }
      else {
	hf_errlog(24040904, "S: EFT reaction: grid format not supported");		
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

///////////////////////////////////////////////////////
#ifdef WITH_PINEAPPL
void RawVec::convolute_PineAPPL() {
  // 1. read the grids if necessary
  // todo: follow Toni's code to deal with exceptions
  if ( ! save_grid_in_memory ) {
    for (string grid_file_name: grid_file_list) {
      pineappl_grid* g = pineappl_grid_read(grid_file_name.c_str());
      pgrid_list.push_back(g);
      hf_errlog(23061202, "I: read PineAPPL grid from " + grid_file_name);
    }
  }
  ////////////////////////////////////////////
  // 2. convolute
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

  // 3. free the grids if necessary
  if (! save_grid_in_memory) {
    for (auto p: pgrid_list)
      pineappl_grid_delete(p);
    pgrid_list.clear();
  }
}
#endif
///////////////////////////////////////////////////////
void RawVec::convolute_APPLgrid() {
  // 1. read the grids if necessary
  if ( ! save_grid_in_memory ) {
    for (string grid_file_name: grid_file_list) {    
	appl::grid* g = new appl::grid(grid_file_name);
	g->trim(); // XMS: why we need this?
	p_APPLgrid_list.push_back(g);
	hf_errlog(24040902, "I: read APPLgrid from " + grid_file_name);
    }
  }
  // 2. convolute
  int shift_bins = 0;

  td->actualizeWrappers(); // XMS: do we need this?
  for (auto pgrid: p_APPLgrid_list) {
    // convolute with all the APPLgrid grids, 
    // and save the results in value_list
    std::vector<double> result = pgrid->vconvolute(
				   pdf_xfxq_wrapper_,
				   pdf_xfxq_wrapper1_,
				   alphas_wrapper_);
    // alphas_wrapper_,
    // order-1,muR,muF,eScale);
    if (shift_bins + result.size() <= value_list.size()) {
      for (size_t i = 0; i < result.size(); ++i)
	value_list[shift_bins + i] = result[i];
    }
    else {
      // error message
    }

  }
  // 3. free the grids if necessary
  if (! save_grid_in_memory) {
  }
}
///////////////////////////////////////////////////////
void RawVec::convolute() {
  /*
    convolute the grids with PDFs.
    if save_grid_in_memory == False, then the grids have to be read into memory
    before every convolution and freed after the convolution.

    results are stored in `value_list`
   */
  if (format == "PineAPPL") {
#ifdef WITH_PINEAPPL
    convolute_PineAPPL();
#endif
  }
  else if (format == "APPLgrid") {
    hf_errlog(24040904, "S: EFT.convolute: grid format not support.");
    convolute_APPLgrid();
  }
  else {
    hf_errlog(24040904, "S: EFT.convolute: grid format not support.");
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
