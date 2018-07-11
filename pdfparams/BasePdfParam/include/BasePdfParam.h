#pragma once


#include <string>


/**
  @class BasePdfParam

  @brief Base clase to deveop PDF parameterisations

  Contains methods to compute PDF value at given x as well as integrals

  @version 0.1
  @date 2018-07-11
  */

class BasePdfParam {
 public:

  /// Default constructor. Name is the PDF name
  BasePdfParam(const std::string& inName){name = inName;}

  /// Compute xf(x,pars). Pure virtual method
  virtual double compute( double const x, double const* pars) = 0;

  /// Compute moments of xf(x) ( 0 -- > of f(x) ), needed for sun-rules
  virtual double moment( double const* pars, int const iMoment = 1) ;
 private:

  /// Name of the PDF object (e.g. uv, dv ...)
  std::string name;

  /// Vector of parameters. Parameters are pointers to doubles
};
