#pragma once

#include <string>
#include <vector>

/**
  @class BaseMinimizer

  @brief Base class to develop interfaces to mimimizers.

  Contains methods to arrange mimimization loop (including lhapdf profiling)

  @version 0.1
  @date 2018-07-25
  */

namespace xfitter
{
  const unsigned int isFreePar  = 1;
  const unsigned int isFixedPar = 2;
  const unsigned int isBoundOar = 3;
    

  
  class BaseMinimizer {
  public:
    /// default constructor
    BaseMinimizer(const std::string& name) :        _name(name),
      _allParameters(),
      _allParameterNames(),
      _fitParameters()
	{}

    /// Initialization
    virtual void initAtStart() = 0;

    /// Provide some information
    virtual void printInfo(){};

    /// Add parameter blocks with Npar parameters. Optional flags and bounds can be used for fixed/bounded parameters.
    virtual void addParameterBlock(int Npar, double const* pars, std::string const* names,  unsigned int const* flags = nullptr, double const* const* bounds = nullptr );

    /// Add single parameter. Optional flags and bounds can be used for fixed/bounded parameters.
    virtual void addParameter(double par, std::string const &name,  unsigned int flag = 0, double  const* bounds = nullptr );
    
    /// Action at each iteration
    virtual void initAtIteration(){};

    /// Miniminzation loop
    virtual void doMimimization() = 0;
    
    /// Perform optional action when equivalent to minuit fcn3 is called (normally after fit)
    virtual void actionAtFCN3(){};

    /// Perform post-fit error analysis
    virtual void errorAnalysis(){};

    /// Number of parameters:
    unsigned int getNAllpars() const { return _allParameters.size(); }

    /// Number of fitted parameters:
    unsigned int getNFitpars() const { return _fitParameters.size(); }
    
  protected:
    /// name to ID minimizer
    std::string _name;

    /// list of all parameters
    std::vector<double> _allParameters;

    /// names of the parameters
    std::vector<std::string> _allParameterNames;
    
    /// list of minimized parameters. Points to sub-set of all parameters
    std::vector<double*> _fitParameters;

  };

  /// For dynamic loader
  typedef BaseMinimizer* create_minimizer();

}
