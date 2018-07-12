#pragma once

#include <string>
#include <map>
#include <functional>

/**
  @class BasePdfDecomposition

  @brief Base class to develop PDF decompositions

  Contains methods to tranform from parameterisation to phyiscal basis

  @version 0.1
  @date 2018-07-11
  */

namespace xfitter
{
  class BasePdfDecomposition {
  public:

    /// Default constructor. Name is the PDF name
    BasePdfDecomposition(const std::string& inName): _name(inName) { };

    /// Optional initialization at the first call
    virtual void initAtStart(const std::string& pars) const = 0;
  
    /// Compute PDF in a physical base in LHAPDF format for given x and Q
    virtual std::function<std::map<int,double>(const double& x)> f0() const = 0;

  private:
    /// Name of the PDF object (e.g. uv, dv ...)
    std::string _name;

  /// Vector of parameters. Parameters are pointers to doubles
  };
}
