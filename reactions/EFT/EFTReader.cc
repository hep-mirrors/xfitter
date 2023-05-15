#include <string>
#include <cstring>
#include "EFTReader.h"
#include "ReactionTheory.h"
#include <yaml-cpp/yaml.h>
#include "EFTReader.h"

using namespace std;

//------------------------------------------------------------------------------------

void EFTReader::setinit(vector<string> name_EFT_param)
{
  num_param = name_EFT_param.size();
  num_bin = 0;

  if (num_param > MAX_NUM_PARAM)
    hf_errlog(23040302, "E: too many EFT parameters");

  // read coefficients for all files
  for (string fname : filename_list){

    YAML::Node coeff_node = YAML::LoadFile(fname);
    int num_bin_one_file = -1;

    // read linear coefficients
    for (int i=0; i < num_param; i++){
      string param_name = name_EFT_param[i];
      if ( coeff_node[param_name] ){
	if (num_bin_one_file < 0) {
	  num_bin_one_file = coeff_node[param_name].size();
	} else if (coeff_node[param_name].size() != num_bin) {
	  hf_errlog(23032903, "E: number of cefficients is not equal to number of bins");
	}
	hf_errlog(23040603, "I: linear coefficients found for: " + param_name);
	num_bin += num_bin_one_file;
	if (coeff[i+1]) {
	  for (double val: coeff_node[param_name].as<std::vector<double> >() ) {
	    (*coeff[i+1]).push_back(val);
	  }
	} else {
	  coeff.insert(std::make_pair(i+1, new vector<double>(coeff_node[param_name].as<std::vector<double> >() )));	  
	}
      } else {
	hf_errlog(23032901, "I: EFT coefficients missing for: " + param_name);
      }
    }

    // read quadratic coefficients
    for (int i=0; i < num_param; i++)
      for (int j=i; j < num_param; j++){
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
	  if (coeff[(i+1)*100 + j+1]) {
	    for (double val : (*pvd) ) (*coeff[(i+1)*100 + j+1]).push_back(val);
	  } else {
	    coeff.insert(std::make_pair((i+1)*100+j+1, pvd));
	  }
	}
      }
  }
}

vector<double> EFTReader::calcxsec(void)
{
  // double xsec[100] = {0.0}; // initialize all xsec as zeros
  // a7: how to save the time of memory allocation? use the reserve method?
  // TODO: parallel calc. for vector operation ?
  //  if (num_bin > 100)
  // hf_errlog(23033001, "E: too many bins");
  vector<double> xsec;
  for (int k=0; k < num_bin; k++)
    xsec.push_back(1.0);

  for (int i=0; i < num_param; i++)
    for (int k=0; k < num_bin; k++)
      if ( coeff[i+1] ){
	xsec[k] += (*coeff[i+1])[k] * val_EFT_param[i];
      }
  
  for (int i=0; i < num_param; i++) 
    for (int j=i; j < num_param; j++)
      for (int k=0; k < num_bin; k++)
	if ( coeff[(i+1)*100+j+1] ) {
	  xsec[k] += (*coeff[(i+1)*100+j+1])[k] * val_EFT_param[i] * val_EFT_param[j]; // / power(labmda, 4);
	}

  // vector<double> vec_xsec;
  // for (int k=0; k < num_bin; k++)
  //   vec_xsec.push_back(xsec[i]);
  // return vec_xsec;
  return xsec;
}

