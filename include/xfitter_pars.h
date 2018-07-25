#pragma once

#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <functional>


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

namespace xfitter{
  // to be defined in evolutions/
  class  BaseEvolution;
  // to be defined in pdfdecompositions/
  class  BasePdfDecomposition;
}

namespace XFITTER_PARS {

  /// Global map of double parameters. They can be used by the minimizer. Initialized based on parameters.yaml
  extern map<string,double*> gParameters;

  /// Global map of integer parameters. Initialized based on parameters.yaml.
  extern map<string,int>    gParametersI;

  /// Global map of string parameters. Initialized based on parameters.yaml.
  extern map<string,string> gParametersS;

  /// Global map of vector parameters. Initialized based on parameters.yaml.
  extern map<string,vector<double> > gParametersV;

  /// Global map of Yaml nodes parameters.. Initialized based on parameters.yaml.
  extern map<string,YAML::Node> gParametersY;

  /// Global map of PDF functions produced by evolutions. 
  extern map<string,std::function<void(double const& x, double const& Q, double* pdfs)> > gXfxQArrays;

  /// Global map to store evolutions
  extern map<string,xfitter::BaseEvolution*> gEvolutions;

  /// Global map to store evolutions
  extern map<string,xfitter::BasePdfDecomposition*> gPdfDecompositions;

  /// Helper function to get string parameters
  std::string getParameterS(std::string name);
  
  /// Parse yaml file @param name
  void parse_file(const std::string& name);

  /// Parse @param node and return maps
  void parse_node(const YAML::Node& node, 
		  std::map<string,double*>& dMap, 
		  std::map<string,int>& iMap, 
		  std::map<string,string>& sMap, 
		  std::map<string,vector<double> >& vMap,
		  std::map<string,YAML::Node> & yMap );


  /// Helper function to register PDFs on the map
  void registerXfxQArray(const string& name, std::function<void(double const& x, double const& Q, double* pdfs)>  xfxArray);

  /// Helper function to retrieve PDFs from the map
  const std::function<void(double const& x, double const& Q, double* pdfs)>  retrieveXfxQArray(const std::string& name);
}

