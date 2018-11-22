
#pragma once

#include "BaseMinimizer.h"

/**
  @class CERESMinimizer

  @brief A class for CERES

  @version 0.1
  @date 2018-08-17
  */



namespace xfitter {
  class CERESMinimizer : public BaseMinimizer
  {
  public:
    /// Default constructor.
    CERESMinimizer () ;

    /// Default constructor. Name is the PDF name
    CERESMinimizer (const std::string& inName);

    /// Optional initialization at the first call
    virtual void atStart() override final;

    /// Minimization loop
    virtual void doMinimization() override final;

    /// Action at last iteration
    virtual void actionAtFCN3() override final;

    /// Error analysis
    virtual void errorAnalysis() override final;

    virtual ConvergenceStatus convergenceStatus()override final;
    /// Parameter transfer
    virtual void addParameter(double par, std::string const &name, double step = 0.01, double const* bounds = nullptr , double  const* priors  = nullptr ) override final;
  private:
    ConvergenceStatus convergence_status=ConvergenceStatus::NORUN;
  };
}
