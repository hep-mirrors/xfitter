#include <string>
#include <cstring>
#include "EFTReader.h"
#include "ReactionTheory.h"
#include <yaml-cpp/yaml.h>
#include "EFTReader.h"

using namespace std;

//------------------------------------------------------------------------------------

void EFTReader::setinit(int num_bin_in, vector<string> name_EFT_param)
{
  num_bin = num_bin_in;
  num_param = name_EFT_param.size();
  if (num_param > MAX_NUM_PARAM)
    hf_errlog(23040302, "E: too many EFT parameters");

  YAML::Node coeff_node = YAML::LoadFile(file_EFT);
  // read linear coefficients
  for (int i=0; i < num_param; i++){
    string param_name = name_EFT_param[i];
    if ( coeff_node[param_name] ){
      // todo: check the size the vector, should = num_bin
      if (coeff_node[param_name].size() != num_bin)
	hf_errlog(23032903, "E: number of cefficients is not equal to number of bins");
      else
	coeff.insert(std::make_pair(-1*(i+1), new vector<double>(	\
				  coeff_node[param_name].as<std::vector<double> >() )));
    }
    else
      hf_errlog(23032901, "I: EFT coefficients missing for: " + param_name);
  }
  // read quadratic coefficients
  for (int i=0; i < num_param; i++)
    for (int j=i; j < num_param; j++){
      string param_name1 = name_EFT_param[i] + "*" + name_EFT_param[j];
      string param_name2 = name_EFT_param[j] + "*" + name_EFT_param[i];
      if (coeff_node[param_name1]) {
  	coeff.insert(std::make_pair(i*100+j, new vector<double>(coeff_node[param_name1].as<std::vector<double> >() )));
      } else if (coeff_node[param_name2]) {
  	coeff.insert(std::make_pair(i*100+j, new vector<double>(coeff_node[param_name2].as<std::vector<double> >() )));
      } else {
  	hf_errlog(23032901, "I: EFT coefficients missing for: " + param_name1);
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
    xsec.push_back(0.0);

  for (int i=0; i < num_param; i++) //a7: calculate the linear contribution
    for (int k=0; k < num_bin; k++) //a7: loop over all bins
      if (coeff[-1*(i+1)]){
	xsec[k] += (*coeff[-1*(i+1)])[k] * val_EFT_param[i]; // / power(lambda, 2);
      }
  
  for (int i=0; i < num_param; i++) 
    for (int j=i; j < num_param; j++)  //a7: calculate the quadratic contribution
      for (int k=0; k < num_bin; k++)  //a7: loop over all bins
	if (coeff[i*100+j]){
	  xsec[k] += (*coeff[i*100+j])[k] * val_EFT_param[i] * val_EFT_param[j]; // / power(labmda, 4);
	}

  // vector<double> vec_xsec;
  // for (int k=0; k < num_bin; k++)
  //   vec_xsec.push_back(xsec[i]);
  // return vec_xsec;
  return xsec;
}

