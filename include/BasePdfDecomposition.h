#pragma once

#include <string>
#include <map>

/**
  @class BasePdfDecomposition

  @brief Base class to develop PDF decompositions

  Contains methods to tranform from parameterisation to physical basis

  @version 0.2
  @date 2019-04-22
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
    /// Returns PDFs in the flavor basis for given x
    /*Parton codes are LHAPDF convention:
     *      i  -6 -5 -4 -3 -2 -1 21 1  2  3  4  5  6  22 11 13 15 -11 -13 -15
     * pdfs[i] tb bb cb sb ub db gl d  u  s  c  b  t  ga e  mu tau ae amu atau
     */
    virtual std::map<int,double>xfxMap(double x)const=0;
    /// Get class name. this is used to verify that the correct concrete class is being used
    virtual const char*getClassName()const=0;
  };
}
