   /*
     @file Vec.cc
     @date 2023-03
     @author X.M. Shen <xmshen137@gmail.com>
   */

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

#ifdef WITH_APPLGRID
#include "appl_grid/appl_grid.h"
#endif

#include "TermData.h" // needed by PDF wrappers

// from Toni's script:
#ifdef WITH_PINEAPPL
#include "pineappl_capi.h"
#endif
#include <sstream>


//--------------------------------------------------------------
using namespace std;
/////////////////////////////////////////////////////////////////////////////
// class Vec()
/////////////////////////////////////////////////////////////////////////////
Vec::Vec(int type_in) {
  type = type_in;
  if (type <= 0) {
    hf_errlog(23051602, "F: EFT: dev: unknown Vec type");
  }
}

void Vec::addIng(RawVec* prvec, double coeff) {
  // ingredient* ping = new ingredient { prvec, coeff };
  bool findQ = false;

  if (type == typem) {
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

// void Vec::bookMonomial(valarray<double>& val) {
void Vec::bookMonomial(double val[99]) {
  if (type != typeMonomial) {
    hf_errlog(24071108, "F: EFT: bug found. please report it to the developer");
  };

  double coeff = 1.0;
  for (const auto& p : param_id_power_list) {
    coeff *= pow(val[p.first], p.second);
  }
  book(coeff);
}

void Vec::addParamPower(int param_id, int power) {
  param_id_power_list.push_back(std::make_pair(param_id, power));
};

/////////////////////////////////////////////////////////////////////////////
// class RawVec()
/////////////////////////////////////////////////////////////////////////////
RawVec::RawVec (YAML::Node node, string key, size_t &num_bin_term, string grid_dir, 
          double xi_ren_in, double xi_fac_in, bool save_grid_Q) {
  // key: tag for the current entry; only used for issuing errors
  // num_bin_term: the number of bins of the EFT term. if = 0, then it 
  //   is determined here, otherwise, it serves as a check.

  xi_ren = xi_ren_in;
  xi_fac = xi_fac_in;
  save_grid_in_memory = save_grid_Q;
  entry = key;
  num_bin = num_bin_term;
  ///////////////////////////////////////////////////////
  // read type and format
  if (node["type"]) {
    string typeS =  node["type"].as<string>();
    if      (typeS == "C")
      type = typeC;
    else if (typeS == "l")
      type = typel;
    else if (typeS == "L")
      type = typeL;
    else if (typeS == "q")
      type = typeq;
    else if (typeS == "Q")
      type = typeQ;
    else if (typeS == "m")
      type = typem;
    else if (typeS == "M")
      type = typeM;
    else if (typeS == "Monomial" || typeS == "monomial")
      type = typeMonomial;
    else
      hf_errlog(23061505, "F: invalid type for entry " + key);      
  }
  else {
    hf_errlog(23061505, "F: type not given for entry " + key);
  }

  if (node["format"]) {
    format = node["format"].as<string>();
    // check if the format is valid
    if (format == FAstr) { // use case/select
    }
    else if (format == FRstr) {
      if (type == typeC) 
        hf_errlog(24071001, "W: EFT: central values should not be ratios");
    }
    else if (format == PineAPPL) {
      if (type == typeMonomial)
        hf_errlog(24071003, "W: EFT: monomial input by grids? Please check " + key);
    } 
    else if (format != APPLgrid) {
        hf_errlog(24071002, "W: EFT: support for APPLgrid is not fully tested.");
        std::cout << "EFT reaction: Warning! "
                  << "Support for APPLgrid is not fully tested. "
                  << "We suggest convert APPLgrid into PineAPPL grids using PineAPPL."
                  << std::endl;
    }
    else 
        hf_errlog(24071004, "F: format " + format + " not supported:" + key);
  }
  else
    hf_errlog(23061504, "F: format not given for entry " + key);

  ///////////////////////////////////////////////////////
  // read optional inputs: xiF, xiR, pdg_id
  if (node["xi_ren"]) 
    setScaleRen(node["xi_ren"].as<double>());
  if (node["xi_fac"]) 
    setScaleFac(node["xi_fac"].as<double>());
  if (node["pdg_id"]) 
    setPDGId(node["pdg_id"].as<int>()); 
    // todo: proton(2212) by default. not realized yet.
    // not sure if pid is supported by APPLgrid?

  ///////////////////////////////////////////////////////
  // read xsec (numbers, or path to grids)
  if (node["xsec"]) {
  }
  else {
    hf_errlog(23061502, "F: EFT: xsec(grids) not given for entry " + key);
  } 
  // read cross section
  if (format == FRstr) {
    ratio_list = node["xsec"].as<vector<double> >();

    // determine or check number of bins
    if (num_bin <= 0) {
      num_bin_term = ratio_list.size();
      num_bin = num_bin_term;
      std::cout << "EFT reaction: num_bin = " << num_bin << std::endl;
    }
    else if (ratio_list.size() != num_bin) {
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cout << "Error: EFT:" << ratio_list.size() << " v.s. " <<  num_bin  << endl;
      hf_errlog(23061501, "F: number of bins does not match for entry " + key);
    }

    // https://cplusplus.com/reference/vector/vector/resize/
    value_list.resize(num_bin, 0.0);
  } // end of FR
  else if (format == FAstr) {
    value_list = node["xsec"].as<vector<double> >();

    if (num_bin <= 0) {
      num_bin_term = value_list.size();
      num_bin = num_bin_term;
      std::cout << "EFT reaction: num_bin = " << num_bin << std::endl;
    }
    else if (value_list.size() != num_bin) {
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cout << "Error: EFT:" << value_list.size() << " v.s. " <<  num_bin  << endl;
      hf_errlog(23061521, "F: number of bins does not match for entry " + key);
    }
  } // end of FA
  else if (format == PineAPPL  || format == APPLgrid) {
    // read path to grids
    for (string filename: node["xsec"].as<vector<string> >()) {
      std::string filepath = grid_dir + "/" + filename;
      std::ifstream file_stream(filepath);
      if (file_stream.good()) {
        grid_file_list.push_back(filepath);
      } 
      else {
        std::cout << "grid file does not exist:"  << std::endl
            << filepath << std::endl
            << grid_dir << std::endl
            << filename << std::endl
            << std::endl;
        hf_errlog(23091201, "F: grid file does not exist:" + filepath);
      }
    } // read grid filename

    // check/determine num_bin,
    // load all the grids if save_grid_in_memory is true
    int num_bin_grids = 0;
    for (string grid_file_name: grid_file_list) {
      if (format == APPLgrid) {
#ifdef WITH_APPLGRID
        // todo: not tested. Do we need TH1D and reference?
        appl::grid* g = new appl::grid(grid_file_name);
        g->trim();
        num_bin_grids += g->Nobs();

        if (save_grid_in_memory) {
          p_APPLgrid_list.push_back(g);
          hf_errlog(24040902, "I: read APPLgrid from " + grid_file_name);
          hf_errlog(24041104, "I: EFT: support for APPLgrid not fully tested");
        }
        else {
          // todo: free the memory?
        }
#else
	      hf_errlog(24040903, "F: APPLgrid support not available");	
#endif        
      }
      else if (format == PineAPPL) {
#ifdef WITH_PINEAPPL
        pineappl_grid* g = pineappl_grid_read(grid_file_name.c_str());
        num_bin_grids += pineappl_grid_bin_count(g);

        if (save_grid_in_memory) {
          pgrid_list.push_back(g);
          hf_errlog(23061202, "I: read PineAPPL grid from " + grid_file_name);
        }
        else {
          pineappl_grid_delete(g);
        }
#else
	      hf_errlog(24040901, "F: PineAPPL support not available");	
#endif
    	} // pineappl
      else {
        hf_errlog(24040904, "F: EFT: grid of type " + format + " not supported");
      }
    } // loop over grid_file_name

    if (num_bin <= 0) {
      num_bin_term = num_bin_grids;
      num_bin = num_bin_term;
      std::cout << "EFT reaction: num_bin = " << num_bin << std::endl;
    }
    else if (num_bin_grids != num_bin) {
      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      hf_errlog(24041501, "F: number of bins does not match for entry " + key);
    }

    value_list.resize(num_bin, 0.0);
    } // end of PineAPPL/APPLgrid
  else {
    hf_errlog(23061503, "F: grid format not support for entry " + key);
  }
  // end of read xsec/grids
  ///////////////////////////////////////////////////////
  if (type != typeC)    
    initReadEFTParam(node);
} // end of constructor

void RawVec::initReadEFTParam(YAML::Node node){
  // check if `param` is provided
  if (type != typeC) {
    if (node["param"]) {
    }
    else
      hf_errlog(23061511, "F: EFT: `param` not found for entry:" +entry);
  }
  else {
    hf_errlog(24071007, "F: assumed to be visible for developers only");
  }

  if (type == typem || type == typeM) {
    // names of parameters
    vector<string> params = node["param"].as<vector<string> >();
    if (params.size() != 2) 
      hf_errlog(23061510, "F: EFT: check `param` for entry:" +entry);
    else {
      param_name1 = params[0];
      param_name2 = params[1];
    }
    // value of parameters
    if (type == typeM) {
      if (node["param_value"]) {
        vector<double> param_vals = node["param_value"].as<vector<double> >();

        if (param_vals.size() != 2) 
          hf_errlog(23061513, "F: `param_value` for entry:" +entry);
        else {
          param_val1 = param_vals[0];
          param_val2 = param_vals[1];
        }
      } 
      else
      	hf_errlog(23061511, "F: `param_value` not found for entry:" +entry);
    }
  } // end of m, M
  else if (type == typeMonomial) {
    // read param names
    param_name_list = node["param"].as<vector<string> >();
    // read powers
    if (node["power"]) {
      power_list = node["power"].as<vector<int> >();
    }
    else
      hf_errlog(24070701, "F: `power` not found for entry:" +entry);

    // check if the powers are valid
    int total_power = 0;
    for (auto power: power_list) {
      total_power += power;
      if (power <= 0) {
        std::cout << "EFT reaction: Warning! non-positive power found in an monomial entry: " +entry
                  << std::endl;
      }
    }
    if (total_power <= 2) {
      std::cout << "EFT reaction: Warning! linear/quadratic entries should have type `l` or `q`, not `monomial`!" +entry
                << std::endl;
    }

    if (power_list.size() != param_name_list.size()) {
      hf_errlog(24070702, "F: sizes of `power` and `param` do not match:" +entry);
    }
  }
  else { // l,q,L,Q
    // names of parameter
    param_name1 = node["param"].as<string>();
    // value of parameter
    if (type < 0) {
      if (node["param_value"])
      	param_val1 = node["param_value"].as<double>();
      else
	      hf_errlog(23061511, "F: `param_value` not found for entry:" +entry);
    }
  }
}


/////////////////////////////////

void RawVec::FR2FA(vector<double> val_list_C) {
  // ratio -> absolute value
  if (format == FRstr) {

    if (value_list.size() != val_list_C.size()) 
      hf_errlog(23061515, "F: size does not match");
    else {
      for (size_t i=0; i<val_list_C.size(); i++)
        value_list[i] = ratio_list[i] * val_list_C[i];
    }
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
#ifdef WITH_APPLGRID
void RawVec::convolute_APPLgrid() {
  // 1. read the grids if necessary
  if ( ! save_grid_in_memory ) {
    for (string grid_file_name: grid_file_list) {    
	appl::grid* g = new appl::grid(grid_file_name);
	g->trim(); // todo XMS: Do and why we need this?
	p_APPLgrid_list.push_back(g);
	hf_errlog(24040902, "I: read APPLgrid from " + grid_file_name);
    }
  }
  // 2. convolute
  int shift_bins = 0;

  for (auto pgrid: p_APPLgrid_list) {
    // convolute with all the APPLgrid grids, 
    // and save the results in value_list
    // XMS: todo: add a new order parameter?
    std::vector<double> result = pgrid->vconvolute(
				   pdf_xfxq_wrapper_,
				   pdf_xfxq_wrapper1_,
				   alphas_wrapper_,
				   1,xi_ren,xi_fac,1.0); // order-1,xi_ren,xi_fac,eScale);
    /* 
       XMS: about QCD order, here 1=NLO(with LO)
       ref:  xfitter/deps/applgrid-1.6.32/src/appl_grid.cxx
             xfitter/src/xfitter_cpp_base.cc:OrderMap()

       order=-1 in xfitter by default; but for APPLgrid -1=corrections_at_NLO(w/o LO)
       otherwise LO=1 in xfitter, while for APPLgrid LO=0
    */
    if (shift_bins + result.size() <= value_list.size()) {
      for (size_t i = 0; i < result.size(); ++i)
	value_list[shift_bins + i] = result[i];
      shift_bins += pgrid->Nobs();
    }
    else {
      hf_errlog(24041101, "F: EFT: number of bins larger than expected.");
    }

  }

  // 3. free the grids if necessary
  if (! save_grid_in_memory) {
    // todo
  }
}
#endif

///////////////////////////////////////////////////////
void RawVec::convolute() {
  /*
    convolute the grids with PDFs.
    if save_grid_in_memory == False, then the grids have to be read into memory
    before every convolution and freed after the convolution.

    results are stored in `value_list`
   */
  if (format == PineAPPL) {
#ifdef WITH_PINEAPPL
    convolute_PineAPPL();
#endif
  }
  else if (format == APPLgrid) {
#ifdef WITH_APPLGRID
    convolute_APPLgrid();
#endif
  }
  else {
    hf_errlog(24040904, "F: EFT.convolute: grid format not support.");
  }
}

/////////////////////////////////
void RawVec::increaseXSecInPlace(valarray<double>& xsec) {
  if (xsec.size() != num_bin)
    hf_errlog(23061516, "F: size does not match");
  else {
    for (size_t i=0; i<value_list.size(); i++) {
      xsec[i] += value_list[i] * coeff;
    }
  }
}
/////////////////////////////////
