#pragma once


#include <string>
#include <map>

/**
  @class BasePdfDecomposition

  @brief Base class to develop PDF decompositions

  Contains methods to tranform from parameterisation to phyiscal basis

  @version 0.1
  @date 2018-07-11
  */

class BasePdfDecomposition {
 public:

  /// Default constructor. Name is the PDF name
  BasePdfDecomposition(const std::string& inName);

  ~BasePdfDecomposition();

  /// Optional initialization at the first call
  virtual void initAtStart(const std::string & pars) {};
  
  /// Compute PDF in a physical base in LHAPDF format for given x and Q
  virtual std::map<int,double> compute(const double& x, const double& Q) const;

 private:

  /// Name of the PDF object (e.g. uv, dv ...)
  std::string name;

  /// Vector of parameters. Parameters are pointers to doubles
};
