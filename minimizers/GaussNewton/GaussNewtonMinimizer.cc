/*
   @file GaussNewtonMinimizer.cc
   @date 2026-10-22

   xFitter glue of the GaussNewton minimiser. See README.md.
*/

#include "GaussNewtonMinimizer.h"
#include "bridge/FortranBridge.h"
#include "bridge/TheoryEvaluator.h"
#include "core/PointModel.h"
#include "xfitter_cpp_base.h"

#include <cmath>
#include <cstdio>

extern "C" {
  void iofilenamesmini_();
  void reset_extra_parameters_();
}

namespace xfitter {

/// the class factory, for dynamic loading
extern "C" GaussNewtonMinimizer* create() {
  return new GaussNewtonMinimizer();
}

GaussNewtonMinimizer::GaussNewtonMinimizer() : BaseMinimizer("GaussNewton") {
}

void GaussNewtonMinimizer::atStart() {
  iofilenamesmini_();
  reset_extra_parameters_();  // the MINUIT extra-parameter priors are not used
}

void GaussNewtonMinimizer::doMinimization() {
  gn::FortranBridge bridge;
  bridge.checkSupported();
  gn::TheoryEvaluator evaluator(*this, bridge);

  Eigen::VectorXd p(getNpars());
  double** pars = getPars();
  for (int i = 0; i < p.size(); ++i) p(i) = *pars[i];
  delete[] pars;

  const double chi2 = evaluator.chi2(p, 1);
  const gn::DataSnapshot snapshot = bridge.snapshot();
  const double moduleChi2 = gn::computePointModel(snapshot).chi2() + evaluator.priorChi2(p);
  const double polarisation = bridge.polarisationChi2();
  std::printf("GaussNewton: %d fitted points, %d nuisance sources, %d parameters\n",
              snapshot.points(), snapshot.sources(), static_cast<int>(p.size()));
  std::printf("GaussNewton: chi2 fcn+priors = %.10g, module = %.10g, polarisation = %g\n", chi2,
              moduleChi2, polarisation);

  // theory-only evaluation at the same point must reproduce THEO (history check)
  const Eigen::VectorXd t = evaluator.theory(p);
  std::printf("GaussNewton: max |t(theory only) - t(fcn)| / sigma = %g\n",
              ((t - snapshot.theory).array() / snapshot.sigma.array()).abs().maxCoeff());

  hf_errlog(2026102210, "F: GaussNewton minimisation is not implemented yet");
}

void GaussNewtonMinimizer::errorAnalysis() {
}

ConvergenceStatus GaussNewtonMinimizer::convergenceStatus() {
  return _convergenceStatus;
}

double GaussNewtonMinimizer::getParameterUncertainty(const std::string& name) const {
  int index = parameterIndex(name);
  if (index < 0 || static_cast<size_t>(index) >= _parameterUncertainties.size()) return std::nan("");
  return _parameterUncertainties[index];
}

}
