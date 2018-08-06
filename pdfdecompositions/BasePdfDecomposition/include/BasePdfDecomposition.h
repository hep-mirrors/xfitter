#pragma once

#include <string>
#include <map>
#include <functional>
#include "BasePdfParam.h"

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
  BasePdfDecomposition(const std::string& inName): _name(inName)  { };
    virtual ~BasePdfDecomposition() {};
    
    /// Optional initialization at the first call
    virtual void initAtStart(const std::string& pars) = 0;
  
    /// Compute PDF in a physical base in LHAPDF format for given x and Q
    virtual std::function<std::map<int,double>(const double& x)> f0() const = 0;

    void addParameterisation(const std::string& pname, BasePdfParam* pdfParam) {
      _pdfParams[pname] = pdfParam;
    }
    
    BasePdfParam* getPdfParam(std::string const& name) const {
      return _pdfParams.at(name);
    }
    
  protected:
    /// PDF parameters
      std::map<std::string,BasePdfParam*> _pdfParams;

    
    
  private:
    /// Name of PDF decomposition 
    std::string _name;

  };

  /// For dynamic loader
  typedef BasePdfDecomposition* create_pdfDecomposition();
}
