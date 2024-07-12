#pragma once

#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <functional>
#include <algorithm>

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
  class BaseEvolution;
  // to be defined in pdfdecompositions/
  class BasePdfDecomposition;
  // to be defined in pdfparams/
  class BasePdfParam;
  // to be defined in minimizers/
  class BaseMinimizer;
}

namespace XFITTER_PARS {
  /// Global pointer to the mimimizer
  extern xfitter::BaseMinimizer* gMinimizer;

  /// Globally available YAML node pointing to root of YAML parameters tree, read from parameters.yaml. Might be modified during runtime
  extern YAML::Node rootNode;
  /// Global map of double parameters. They can be used by the minimizer. Initialized based on parameters.yaml
  extern map<string,double*> gParameters;

  /// Global map of integer parameters. Initialized based on parameters.yaml.
  extern map<string,int>    gParametersI;

  /// Global map of string parameters. Initialized based on parameters.yaml.
  extern map<string,string> gParametersS;

  /// Global map of vector double parameters. Initialized based on parameters.yaml.
  extern map<string,vector<double> > gParametersV;

  /// Global map of vector string parameters. Initialized based on parameters.yaml.
  extern map<string,vector<string> > gParametersVS;

  /// Global map of Yaml nodes parameters.. Initialized based on parameters.yaml.
  extern map<string,YAML::Node> gParametersY;
  /// Global map to store evolutions
  /// Do not access it directly, use xfitter::get_evolution
  extern map<string,xfitter::BaseEvolution*> gEvolutions;
  /// Global map to store decompositions
  /// Do not access it directly, use xfitter::get_pdfDecomposition
  extern map<string,xfitter::BasePdfDecomposition*> gPdfDecompositions;
  /// Global map to store parameterisations
  /// Do not access it directly, use xfitter::getParameterisation
  extern map<string,xfitter::BasePdfParam*>gParameterisations;
  /// Helper function to get input decomposition from a yaml node of evolution
  ///
  /// It finds a "decomposition" subnode in given node, extracts a decomposition name from it, finds this decomposition and returns it
  /// This function will either return a valid function, or issue a fatal error and never return
  /// This is used in evolution initialization
  xfitter::BasePdfDecomposition*getInputDecomposition(const YAML::Node&evolutionNode);
  /// Helper function to get a yaml node corresponding to an evolution, by this evolutions's instance name
  YAML::Node getEvolutionNode(const std::string&name="");
  /// Helper function to get a yaml node corresponding to a decomposition, by this decomposition's instance name
  YAML::Node getDecompositionNode(const std::string&name="");
  /// Helper function to get a yaml node corresponding to a parameterisation, by this parameterisation's instance name
  YAML::Node getParameterisationNode(const std::string&name="");

  /// Helper functions
  string getDefaultEvolutionName();
  string getDefaultDecompositionName();
  /// Functions to get parameters from corresponding maps but with better reporting of errors
  double*getParamD(const string&name);
  int    getParamI(const string&name);
  string getParamS(const string&name);

  /// Functions to get parameters from nodes but with better reporting of errors
  string getParamFromNodeS(const string& name, const YAML::Node& node);

  /// Functions to get evolution parameters from YAML steering
  /// They either return the parameter, or report an error and stop the program
  double getEvolutionParamD(const string& evolutionName,const string& parameterName);
  //TODO: same for string and int, maybe?

  /// Parse @param node and return maps
  void parse_node(const YAML::Node& node,
		  std::map<string,double*>& dMap,
		  std::map<string,int>& iMap,
		  std::map<string,string>& sMap,
		  std::map<string,vector<double> >& vMap,
		  std::map<string,vector<string> >& vsMap,
		  std::map<string,YAML::Node> & yMap );
  /// Allocate memory for a new double*-typed parameter, record it in gParameters and return it
  /// This kind of parameter is not passed to minimizer
  const double*createConstantParameter(const string&name,double value);
}

