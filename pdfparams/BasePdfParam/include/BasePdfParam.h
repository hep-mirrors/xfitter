#pragma once


#include <string>
#include <yaml-cpp/yaml.h>

/**
  @class BasePdfParam

  @brief Base class to develop PDF parameterisations

  Contains methods to compute PDF value at given x as well as integrals

  @version 0.1
  @date 2018-07-11
  */

class BasePdfParam {
 public:

  /// Default constructor. Name is the PDF name
  BasePdfParam(const std::string& inName): _name(inName),_npar(0) {}

  /// Set number of parameters
  void SetNPar(int npar) { _npar = npar;}
  
  /// Compute xf(x,pars). Pure virtual method
  virtual double compute( double const x, double const* pars) = 0;

  /// Compute moments of xf(x) ( for i=-1  of f(x) ), needed for sun-rules
  virtual double moment( double const* pars, int const iMoment = 0) ;

  /// Return number of parameters
  const int getNPar() const {return _npar;} 

  /// Get name of the PDF object
  const std::string getName() const {return _name;} 

  /// Get initial values from a yaml node. Uses node[getName] as the basis
  void initFromYaml(const YAML::Node& node, double* pars) ;
  
 private:
  /// Name of the PDF object (e.g. uv, dv ...)
  std::string _name;
  /// Number of parameters
  int _npar;
  
  /// Vector of parameters. Parameters are pointers to doubles
};
