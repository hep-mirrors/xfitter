#pragma once

#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>


using std::map;
using std::string;
using std::vector;

/*!
 @file xfitter_pars.h
 @date Sun 16 April 2017
 @author SG

 Contains functions to read parameters.yaml, 
 global maps to store parameters,  and fortran interface functions.
*/

namespace XFITTER_PARS {

  extern map<string,double*> gParameters;
  extern map<string,int>    gParametersI;
  extern map<string,string> gParametersS;
  extern map<string,vector<double> > gParametersV;
  extern map<string,YAML::Node> gParametersY;


  // Helper functions
  double getParDouble(const string& name);
  
  /// Parse yaml file @param name
  void parse_file(const std::string& name);

  /// Parse @param node and return maps
  void parse_node(const YAML::Node& node, 
		  std::map<string,double*>& dMap, 
		  std::map<string,int>& iMap, 
		  std::map<string,string>& sMap, 
		  std::map<string,vector<double> >& vMap,
		  std::map<string,YAML::Node> & yMap );
  
}

