/*
   @file GaussNewtonMinimizer.cc
   @date 2026-10-22

   xFitter glue of the GaussNewton minimiser. See README.md.
*/

#include "GaussNewtonMinimizer.h"
#include "bridge/FortranBridge.h"
#include "bridge/TheoryEvaluator.h"
#include "core/PointModel.h"
#include "jacobian/FiniteDifferenceJacobian.h"
#include "xfitter_cpp_base.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

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

  // finite-difference Jacobian, in this process and with 2 workers
  const double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd lower = Eigen::VectorXd::Constant(p.size(), -inf);
  Eigen::VectorXd upper = Eigen::VectorXd::Constant(p.size(), inf);
  for (int i = 0; i < p.size(); ++i)
    if (_bounds[i]) {
      lower(i) = _bounds[i][0];
      upper(i) = _bounds[i][1];
    }
  Eigen::VectorXd steps(p.size());
  for (int i = 0; i < p.size(); ++i) {
    const auto it = _steps.find(_allParameterNames[i]);
    steps(i) = it != _steps.end() ? it->second : 0.01 * std::max(std::abs(p(i)), 1.);
  }
  gn::FDOptions options;
  gn::FiniteDifferenceJacobian serial(evaluator, steps, lower, upper, options);
  const Eigen::MatrixXd J = serial.compute(p, snapshot.theory);
  options.workers = 2;
  gn::FiniteDifferenceJacobian forked(evaluator, steps, lower, upper, options);
  const Eigen::MatrixXd Jf = forked.compute(p, snapshot.theory);
  std::printf("GaussNewton: max |J(serial) - J(2 workers)| = %g, |J| = %g, %d evaluations\n",
              (J - Jf).cwiseAbs().maxCoeff(), J.norm(), serial.evaluations());
  // linear prediction for a step of 10 FD steps in all parameters
  const Eigen::VectorXd dp = 10 * options.relativeStep * p.cwiseAbs().cwiseMax(steps);
  const Eigen::VectorXd moved = evaluator.theory(p + dp);
  std::printf("GaussNewton: linear prediction max |t(p+dp) - t - J dp| / sigma = %g (|t(p+dp) - t| / sigma = %g)\n",
              ((moved - snapshot.theory - J * dp).array() / snapshot.sigma.array()).abs().maxCoeff(),
              ((moved - snapshot.theory).array() / snapshot.sigma.array()).abs().maxCoeff());

  hf_errlog(2026102210, "F: GaussNewton minimisation is not implemented yet");
}

void GaussNewtonMinimizer::errorAnalysis() {
}

ConvergenceStatus GaussNewtonMinimizer::convergenceStatus() {
  return _convergenceStatus;
}

void GaussNewtonMinimizer::addParameter(double par, std::string const& name, double step,
                                        double const* bounds, double const* priors) {
  BaseMinimizer::addParameter(par, name, step, bounds, priors);
  if (step > 0) _steps[name] = step;
}

double GaussNewtonMinimizer::getParameterUncertainty(const std::string& name) const {
  int index = parameterIndex(name);
  if (index < 0 || static_cast<size_t>(index) >= _parameterUncertainties.size()) return std::nan("");
  return _parameterUncertainties[index];
}

}
