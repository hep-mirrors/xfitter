#pragma once

#include "BaseMinimizer.h"

#include <map>
#include <string>
#include <vector>

/**
  @class GaussNewtonMinimizer

  @brief Levenberg-Marquardt minimisation of the profiled chi2 with a theory Jacobian

  The minimiser only proposes steps. Every step is accepted or rejected with the official
  chi2 returned by fcn. The mathematics is described in README.md.

  @date 2026-10-22
  */

namespace xfitter {
  class GaussNewtonMinimizer : public BaseMinimizer
  {
  public:
    GaussNewtonMinimizer();

    virtual void atStart() override final;
    virtual void doMinimization() override final;
    virtual void errorAnalysis() override final;
    virtual ConvergenceStatus convergenceStatus() override final;
    virtual double getParameterUncertainty(const std::string& name) const override final;
    virtual void addParameter(double par, std::string const& name, double step = 0.01,
                              double const* bounds = nullptr, double const* priors = nullptr) override final;

  private:
    ConvergenceStatus _convergenceStatus = ConvergenceStatus::NORUN;
    std::vector<double> _parameterUncertainties;
    /// YAML steps by name. Parameters added from Fortran (AddExternalParam) bypass
    /// addParameter and get a default step.
    std::map<std::string, double> _steps;
  };
}
