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

  assert (filename_list.size() == 1);
  string fname = filename_list[0];
  YAML::Node node = YAML::LoadFile(fname);

  for (YAML::const_iterator entry=node.begin(); entry!=node.end(); ++entry ) {
    if (entry["mask"])
      if (entry["mask"].as<bool>() == true)
	continue;

    string type;
    string param_name;
    if (entry["type"])
      type = entry["type"].as<string>(); // use string instead of char for future extension
    else
      hf_errlog(23052303, "F: type not given for entry " + entry.first);

    // C term
    if (type == "C") {
      assert (prvec_C == nullptr);
      prvec_C = new RawVe(entry);
      assert (prvec_C->format != "FR");
    }
    else if (type=="l" || type=="L" || type=="q" || type=="Q") {
      if (entry["param"]) {
	param_name = entry["param"].as<string>();
	if (find_EFT_param_id.count(param_name) == 0) 
	  continue;
	int id = find_EFT_param_id[param_name];
	if (basis.count[id] == 0) {
	  basis.insert = make_pair(id, new Vec(1));
	}

	RawVec* prvec = new RawVec(entry);
	raw_basis.push_back(prvec);
	basis[i].addIng(prvec, 1.0);
      }
      else {
	hf_errlog(23052401, "F: param name not given for " + entry.first);
      }
    } // end of l, q
    else if (type=="m" || type=="M") {
      if (entry["param"]) {
	vector<string> param_names = entry["param"].as<vector<string> >();
	// todo assert we can find param_val
	vector<double> param_vals = entry["param_val"].as<vector<double> >();
	
	assert (param_names.size() == 2);
	if (find_EFT_param_id.count(param_names[0]) * find_EFT_param_id.count(param_names[1]) == 0 ) 
	  continue;

	int i1 = find_EFT_param_id[param_names[0]];
	int i2 = find_EFT_param_id[param_names[1]];
	double param_val1;
	double param_val2;


	if (i1 > i2) {
	  i1 += i2;
	  i2 = i1 - i2;
	  i1 -= i1;
	  param_val1 = param_vals[1];
	  param_val2 = param_vals[0];	  
	}
	else{
	  param_val1 = param_vals[0];
	  param_val2 = param_vals[1];	  
	}

	assert(basis.count(i1*100+i2) == 0);
	basis.insert(make_pair(i1*100+i2, new Vec(2)));
	
	RawVec* prvec = new RawVec(entry);
	raw_basis.push_back(prvec);
	basis[i1*100+i2].addIng(prvec, 1.0);
      }
    } // end of m, M
    else {
      hf_errlog(23052402, "F: type not supported for entry " + entry.first);
    }
  } // end of loop over entries

  initlq();
  initm();
  initrvec();

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
  if (inputType == "fixed")
    return calcXSecFixed();
  else
    return calcXSecMixed();
};

//------------------------------------------------------------------------------------
vector<double> EFTReader::calcXSecMixed(){
  vector<double> xsec;

  // C
  for (int k=0; k < num_bin; k++)
    xsec.push_back(0.0);
  
  if (no_central == false) {
    prvec_C->increase_xsec_in_place(xsec); // to do, pass the referene instead !!
  }

  // lqm
  for (auto prvec : raw_basis) 
    prvec->increase_xsec_in_place(xsec);

  if (abs_output == false)
    for (int i=0; i<num_bin; ++i)
      xsec[i] /= prvec_C->value_list[i];

  return xsec;
};

//------------------------------------------------------------------------------------
vector<double> EFTReader::calcXSecFixed(){

  if (debug > 0) {
    std::cout << "=======================================================" << std::endl;
    std::cout << "EFTReader.calcXSec: size of map coeff: " << coeff.size() << std::endl;
  }

  vector<double> xsec;

  // C
  if (no_central == true) {
    for (int k=0; k < num_bin; k++)
      xsec.push_back(0.0);
  }
  else {
    for (int k=0; k < num_bin; k++)
      xsec.push_back(1.0);
  }

  // l
  for (int i=0; i < num_param; i++) {
    if (coeff.count(i+1) > 0) {

      for (int k=0; k < num_bin; k++) {
	xsec[k] += (*(coeff[i+1]))[k] * val_EFT_param[i];
      }
      /////////////////////////////////
      if (debug > 0) {
	std::cout << "=======================================================" << std::endl;
	std::cout << "EFTReader.calcXSec: l:" + name_EFT_param[i] + " = " << val_EFT_param[i] << std::endl;
	std::cout << "EFTReader.calcXSec: new xsec=" << std::endl;
	for (double v : xsec) {
	  std::cout << v << ", ";
	}
	std::cout << std::endl;
      }
      /////////////////////////////////
    } 
    else {
      hf_errlog(23051601, "F: linear coefficients missing for: " + name_EFT_param[i]);
    }
  
  }

  // q&m
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
  } // end of q&m

  return xsec;
}
