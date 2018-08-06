#pragma once

#include <string>
#include <map>

/**
  @class BaseMinimizer

  @brief Base class to develop interfaces to mimimizers.

  Contains methods to arrange mimimization loop (including lhapdf profiling)

  @version 0.1
  @date 2018-07-25
  */

namespace xfitter
{
  class BaseMinimizer {
  public:
    /// default constructor
    BaseMinimizer(const std::string& name) :        _name(name),
      _fitParameters(nullptr),
      _allParameters(nullptr) {}

    /// Initialization
    virtual void initAtStart() = 0;

    /// Provide some information
    virtual void printInfo(){};

    /// Add parameter blocks with Npar parameters. Optional flags and bounds can be used for fixed/bounded parameters.
    virtual void addParameterBlock(int Npar, double const* pars, unsigned int const* flags = nullptr, double const* const* bounds = nullptr );
    
    /// Action at each iteration
    virtual void initAtIteration(){};

    /// Miniminzation loop
    virtual void doMimimization() = 0;
    
    /// Perform optional action when equivalent to minuit fcn3 is called (normally after fit)
    virtual void actionAtFCN3(){};

    /// Perform post-fit error analysis
    virtual void errorAnalysis(){};
    
  private:
    /// name to ID minimizer
    std::string _name;

    /// list of minimized parameter blocks
    double** _fitParameters;

    /// list of all parameter blocs
    double** _allParameters;

    /// maping from fitted parameters to all parameters
  };

  /// For dynamic loader
  typedef BaseMinimizer* create_minimizer();

}
