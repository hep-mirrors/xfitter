
#pragma once

#include "BaseMinimizer.h"

/**
  @class MINUITMinimizer

  @brief A class for MINUIT

  @version 0.1
  @date 2018-08-17
  */

namespace xfitter {

class MINUITMinimizer : public BaseMinimizer
{
  public:
     /// Default constructor.
    MINUITMinimizer () ;

     /// Default constructor. Name is the PDF name
    MINUITMinimizer (const std::string& inName);

    /// Optional initialization at the first call
    virtual void atStart() override final;

    /// Minimization loop
    virtual void doMinimization() override final;

    /// Action at last iteration
    virtual void actionAtFCN3() override final;

    /// Error analysis
    virtual void errorAnalysis() override final;

    virtual ConvergenceStatus convergenceStatus()override final;
    /// parameters
    virtual void addParameter(double par, std::string const &name, double step = 0.01, double const* bounds = nullptr , double  const* priors  = nullptr ) override final;
    virtual double getParameterUncertainty(const std::string& name) const override final;
    virtual double getParameterValue(const std::string& name) const override final;
    virtual double getParameterValue(const std::string& name, const double* pars) const override final;
    virtual unsigned int getNfreePars() const override final;
};
}
