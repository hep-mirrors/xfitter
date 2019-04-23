#pragma once

#include <string>
#include <map>
#include <functional>
#include "BasePdfParam.h"
// #include <iostream>

/**
  @class BasePdfDecomposition

  @brief Base class to develop PDF decompositions

  Contains methods to tranform from parameterisation to physical basis

  @version 0.1
  @date 2018-07-11
  */

namespace xfitter
{
  class BasePdfDecomposition {
  public:
		const std::string _name;//Unique name used to identify this decomposition instance
    /// Default constructor.
    BasePdfDecomposition(const char*name):_name(name){}
    virtual ~BasePdfDecomposition(){}

    /// Initialization at the first call
    virtual void atStart(){}
    /// Optional initialization at each iteration. Can be used to compute sum-rules
    virtual void atIteration(){}
    /// This function should be called when at least one parameter in the YAML node of given decomposition changes
    virtual void atConfigurationChange(){}
    /// Returns a LHAPDF-style function, that returns PDFs in a physical basis for given x
    virtual std::function<std::map<int,double>(const double& x)> f0() const = 0;
    /// Get class name, can be used to verify that the correct concrete class is being used
    virtual const char*getClassName()const=0;
  };
  /// For dynamic loader
  //Wait, is this even used anywhere? --Ivan
  typedef BasePdfDecomposition* create_pdfDecomposition();
}
