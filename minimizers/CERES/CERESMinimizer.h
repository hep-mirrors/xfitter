
#pragma once

#include "BaseMinimizer.h"
#include "ceres/ceres.h"

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

    //Output functions
    void writePars(const double* covmat);
    void writeOutput(ceres::Solver::Summary mySummary, const double* covmat);

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

    static double glboff;
    static double *offset;
    static double totoffset;
    //the variables offset and totoffset are not thread safe

  private:
    ConvergenceStatus convergence_status=ConvergenceStatus::NORUN;
  };
}
