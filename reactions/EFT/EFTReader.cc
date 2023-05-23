#include <string>
#include <cstring>
#include "EFTReader.h"
#include <yaml-cpp/yaml.h>
#include "EFTReader.h"

using namespace std;

//------------------------------------------------------------------------------------
void EFTReader::initParamName(vector<string> name_EFT_param_in){
  int i = 1;
  num_param = name_EFT_param_in.size();

  if (num_param > MAX_NUM_PARAM)
    hf_errlog(23040302, "E: too many EFT parameters");

  for (string name : name_EFT_param_in) {
    name_EFT_param.push_back(name);
    find_EFT_param_id.insert(std::make_pair(name, i++));
  }

  assert(i == num_param-1);
}

//------------------------------------------------------------------------------------
void EFTReader::read_input(){
  num_bin = 0;  
  if (inputType == "fixed") 
    read_fixed_input();
  else
    read_mixed_input();
}

//------------------------------------------------------------------------------------
void EFTReader::read_fixed_input() {

  // read coefficients for all files
  for (string fname : filename_list) {

    YAML::Node coeff_node = YAML::LoadFile(fname);
    int num_bin_one_file = -1;

    // read linear coefficients
    for (int i=0; i < num_param; i++) {
      string param_name = name_EFT_param[i];

      if (coeff_node[param_name]) {
	// check the number of bins
	if (num_bin_one_file < 0) {
	  num_bin_one_file = coeff_node[param_name].size();
	  num_bin += num_bin_one_file;
	} 
	else if (coeff_node[param_name].size() != num_bin_one_file) {
	  hf_errlog(23032903, "F: # coefficients != # bins");
	}

	// read the coeff.
	if (coeff.count(i+1) > 0) {
	  for (double val: coeff_node[param_name].as<std::vector<double> >() ) {
	    (*coeff[i+1]).push_back(val);
	  }
	} 
	else {
	  coeff.insert(std::make_pair(i+1, new vector<double>(coeff_node[param_name].as<std::vector<double> >() )));

	  if (debug > 0) {
	    std::cout << "=======================================================" << std::endl;
	    std::cout << "EFTReader.init: size of map coeff: " << coeff.size() << std::endl;
	  }
	}
      } 
      else {
	hf_errlog(23032901, "F: EFT coefficients missing for: " + param_name);
      }
    } // end of reading linear coeff.

    // read quadratic coefficients
    for (int i=0; i < num_param; i++) {
      for (int j=i; j < num_param; j++) {
	string param_name1 = name_EFT_param[i] + "*" + name_EFT_param[j];
	string param_name2 = name_EFT_param[j] + "*" + name_EFT_param[i];
	vector<double>* pvd;
	bool found = false;
	if (coeff_node[param_name1]) {
	  pvd = new vector<double>(coeff_node[param_name1].as<std::vector<double> >());
	  found = true;
	} else if (coeff_node[param_name2]) {
	  pvd = new vector<double>(coeff_node[param_name2].as<std::vector<double> >());
	  found = true;
	} else {
	  hf_errlog(23032901, "I: EFT coefficients missing for: " + param_name1);
	}

	if (found) {
	  if (coeff.count((i+1)*100 + j+1) > 0) {
	    for (double val : (*pvd)) (*coeff[(i+1)*100 + j+1]).push_back(val);
	  } else {
	    coeff.insert(std::make_pair((i+1)*100+j+1, pvd));
	  }
	}
      }
    } // end of reading quadratic coeff.

  }
}

//------------------------------------------------------------------------------------
void EFTReader::read_mixed_input(){

}
//------------------------------------------------------------------------------------

void EFTReader::initlq(){
}
//------------------------------------------------------------------------------------

void EFTReader::initm(){
}

//------------------------------------------------------------------------------------
void EFTReader::initrvec(){
}

//------------------------------------------------------------------------------------
void EFTReader::initIter(vector<double> list_val){
  setValEFT(list_val);
  
  if (inputType == "mixed") {
    updatervec();
    book();
  }
}

//------------------------------------------------------------------------------------
void EFTReader::updatervec(); {

} 

//------------------------------------------------------------------------------------
void EFTReader::book(){

}

//------------------------------------------------------------------------------------
void EFTReader::setValEFT(vector<double> list_val) {
  // executed for each computation
  if (num_param == list_val.size()) {
    for (int i=0; i<num_param; i++)
      val_EFT_param[i] = list_val[i];
  } else {
    hf_errlog(23040301, "E: number of EFT parameters does not match");
  }

  if (debug > 2) {
    std::cout << "=======================================================" << std::endl;
    std::cout << "EFTReader.setValEFT" << std::endl;
    for (int i=0; i<num_param; i++) {
      std::cout << name_EFT_param[i] << "=" <<  val_EFT_param[i] << std::endl;
    }
  }
};

//------------------------------------------------------------------------------------
vector<double> EFTReader::calcXSec() {

};

//------------------------------------------------------------------------------------
vector<double> EFTReader::calcXSecMixed(){
};

//------------------------------------------------------------------------------------
vector<double> EFTReader::calcXSecFixed(){
  // double xsec[100] = {0.0}; // initialize all xsec as zeros
  // a7: how to save the time of memory allocation? use the reserve method?
  if (debug > 0) {
    std::cout << "=======================================================" << std::endl;
    std::cout << "EFTReader.calcXSec: size of map coeff: " << coeff.size() << std::endl;
  }

  vector<double> xsec;
  for (int k=0; k < num_bin; k++)
    xsec.push_back(1.0);

  for (int i=0; i < num_param; i++) {
    if (coeff.count(i+1) > 0) {

      for (int k=0; k < num_bin; k++) {
	xsec[k] += (*(coeff[i+1]))[k] * val_EFT_param[i];
      }

      if (debug > 0) {
	std::cout << "=======================================================" << std::endl;
	std::cout << "EFTReader.calcXSec: l:" + name_EFT_param[i] + " = " << val_EFT_param[i] << std::endl;
	std::cout << "EFTReader.calcXSec: new xsec=" << std::endl;
	for (double v : xsec) {
	  std::cout << v << ", ";
	}
	std::cout << std::endl;
      }

    } else {
      hf_errlog(23051601, "F: linear coefficients missing for: " + name_EFT_param[i]);
    }
  
  }

  for (int i=0; i < num_param; i++) {
    for (int j=i; j < num_param; j++) {
      if (coeff.count((i+1)*100+j+1) > 0) {

	vector<double>* pvec = coeff[(i+1)*100+j+1];
	for (int k=0; k < num_bin; k++)
	  xsec[k] += (*pvec)[k] * val_EFT_param[i] * val_EFT_param[j];
	/////////////////////////////////
	if (debug > 0) {
	  std::cout << "=======================================================" << std::endl;
	  std::cout << "EFTReader.calcXSec: q/m: " << name_EFT_param[i] +  "*" + name_EFT_param[j] 
		    << " = " << val_EFT_param[i] * val_EFT_param[j] << std::endl;
	  std::cout << "EFTReader.calcXSec:" << std::endl;
	  for (double v : xsec) {
	    std::cout << v << ", ";
	  }
	  std::cout << std::endl;
	}
	/////////////////////////////////
      }
    }
  } // end of double loop

  return xsec;
}
