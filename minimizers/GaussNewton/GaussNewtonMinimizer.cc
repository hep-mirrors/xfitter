/*
   @file GaussNewtonMinimizer.cc
   @date 2026-10-22

   xFitter glue of the GaussNewton minimiser: steering, wiring of the pieces, outputs.
   See README.md.
*/

#include "GaussNewtonMinimizer.h"
#include "DerivativeCheck.h"
#include "ProfiledChi2.h"
#include "bridge/FortranBridge.h"
#include "bridge/TheoryEvaluator.h"
#include "core/LevenbergMarquardt.h"
#include "core/ObservedHessian.h"
#include "jacobian/FiniteDifferenceJacobian.h"
#include "util/Report.h"
#include "util/Timer.h"
#include "xfitter_cpp_base.h"
#include "xfitter_pars.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <set>

extern "C" {
  void iofilenamesmini_();
  void reset_extra_parameters_();
  void errbandssym_();
}

namespace xfitter {

/// the class factory, for dynamic loading
extern "C" GaussNewtonMinimizer* create() {
  return new GaussNewtonMinimizer();
}

namespace {

struct Settings {
  gn::LMOptions lm;
  gn::FDOptions fd;
  bool covariance = true;
  bool observedHessian = true;  ///< covariance from second differences of the chi2 (as HESSE)
  bool checkDerivatives = false;
};

/// The GaussNewton block of parameters.yaml (README.md, "Steering")
Settings readSettings() {
  Settings s;
  const YAML::Node node = XFITTER_PARS::rootNode["GaussNewton"];
  if (!node.IsMap()) return s;
  const std::set<std::string> known = {"maxIterations", "tolerance", "lambda0", "jacobian", "fdScheme",
                                       "fdStep", "threads", "covariance", "hessian", "doErrors",
                                       "checkDerivatives"};
  for (const auto& item : node)
    if (!known.count(item.first.as<std::string>()))
      hf_errlog(2026102250, "W: GaussNewton: unknown steering key " + item.first.as<std::string>());

  s.lm.maxIterations = node["maxIterations"].as<int>(s.lm.maxIterations);
  s.lm.tolerance = node["tolerance"].as<double>(s.lm.tolerance);
  s.lm.lambda0 = node["lambda0"].as<double>(s.lm.lambda0);
  if (node["jacobian"].as<std::string>("FiniteDifference") != "FiniteDifference")
    hf_errlog(2026102251, "F: GaussNewton: only jacobian: FiniteDifference is available");
  const std::string scheme = node["fdScheme"].as<std::string>("forward");
  if (scheme != "forward" && scheme != "central")
    hf_errlog(2026102252, "F: GaussNewton: fdScheme must be forward or central");
  s.fd.central = scheme == "central";
  s.fd.relativeStep = node["fdStep"].as<double>(s.fd.relativeStep);
  s.fd.workers = node["threads"].as<int>(s.fd.workers);
  s.covariance = node["covariance"].as<int>(1) != 0;
  const std::string hessian = node["hessian"].as<std::string>("observed");
  if (hessian != "observed" && hessian != "gaussnewton")
    hf_errlog(2026102259, "F: GaussNewton: hessian must be observed or gaussnewton");
  s.observedHessian = hessian == "observed";
  s.checkDerivatives = node["checkDerivatives"].as<int>(0) != 0;
  return s;
}

}  // namespace

GaussNewtonMinimizer::GaussNewtonMinimizer() : BaseMinimizer("GaussNewton") {
}

void GaussNewtonMinimizer::atStart() {
  iofilenamesmini_();
  reset_extra_parameters_();  // the MINUIT extra-parameter priors are not used
}

void GaussNewtonMinimizer::doMinimization() {
  const auto start = std::chrono::steady_clock::now();
  const Settings settings = readSettings();
  const int npar = getNpars();

  gn::FortranBridge bridge;
  bridge.checkSupported();
  gn::TheoryEvaluator evaluator(*this, bridge);

  Eigen::VectorXd p0(npar), steps(npar);
  Eigen::VectorXd lower = Eigen::VectorXd::Constant(npar, -std::numeric_limits<double>::infinity());
  Eigen::VectorXd upper = Eigen::VectorXd::Constant(npar, std::numeric_limits<double>::infinity());
  double** pars = getPars();
  for (int i = 0; i < npar; ++i) {
    p0(i) = *pars[i];
    const auto it = _steps.find(_allParameterNames[i]);
    steps(i) = it != _steps.end() ? it->second : 0.01 * std::max(std::abs(p0(i)), 1.);
    if (_bounds[i]) {
      lower(i) = _bounds[i][0];
      upper(i) = _bounds[i][1];
    }
  }
  delete[] pars;

  gn::FiniteDifferenceJacobian jacobian(evaluator, steps, lower, upper, settings.fd);
  gn::ProfiledChi2 problem(evaluator, bridge, jacobian, lower, upper);

  // First call (fcn flag 1) and start-point checks
  problem.setFcnFlag(1);
  problem.chi2(p0);
  if (bridge.polarisationChi2() != 0)
    hf_errlog(2026102253, "F: GaussNewton: polarisation shifts are not supported");
  const gn::DataSnapshot first = bridge.snapshot();
  const Eigen::VectorXd t = evaluator.theory(p0);
  if (((t - first.theory).array() / first.sigma.array()).abs().maxCoeff() > 1e-6)
    hf_errlog(2026102254, "W: GaussNewton: theory depends on the evaluation history "
                          "(e.g. a dependent parameter using a sum-rule normalisation)");
  std::printf("GaussNewton: %d fitted points, %d nuisance sources, %d parameters\n", first.points(),
              first.sources(), npar);

  if (settings.checkDerivatives) {
    Eigen::VectorXd h(npar);
    for (int i = 0; i < npar; ++i) h(i) = settings.fd.relativeStep * std::max(std::abs(p0(i)), steps(i));
    gn::checkDerivatives(problem, p0, h, _allParameterNames);
  }

  // Minimisation, Eq. (11)
  std::printf("%5s %14s %11s %11s %9s %9s %6s %6s %8s %8s\n", "iter", "chi2", "decrement", "predicted",
              "rho", "lambda", "trials", "fcn", "theory", "time[s]");
  const gn::LMResult result = gn::LevenbergMarquardt(settings.lm).minimise(
      problem, p0, [&](const gn::LMIteration& it) {
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::printf("%5d %14.6f %11.4e %11.4e %9.3f %9.2e %6d %6d %8d %8.1f\n", it.iteration, it.chi2,
                    it.decrement, it.predicted, it.rho, it.lambda, it.trials, evaluator.fcnCalls(),
                    jacobian.evaluations(), elapsed);
        std::printf("      p =");
        for (int i = 0; i < it.x.size(); ++i) std::printf(" %.6g", it.x(i));
        std::printf("\n");
        std::fflush(stdout);
      });
  std::printf("GaussNewton: %s after %d iterations, chi2 = %.6f\n", result.message.c_str(),
              result.iterations, result.chi2);

  // Covariance, Eq. (10)
  Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(npar, npar);
  bool positive = true;
  int hessianCalls = 0;
  if (settings.covariance) {
    Eigen::MatrixXd hessian = result.model.hessian;
    if (settings.observedHessian) {
      const gn::ObservedHessianResult observed = gn::observedHessian(problem, result.x, result.chi2, hessian);
      hessianCalls = observed.chi2Calls;
      Eigen::MatrixXd gnCovariance;
      gn::covarianceFromHessian(hessian, gnCovariance);
      std::printf("GaussNewton: observed Hessian from %d chi2 calls; errors observed / Gauss-Newton:", hessianCalls);
      hessian = observed.hessian;
      Eigen::MatrixXd observedCovariance;
      gn::covarianceFromHessian(hessian, observedCovariance);
      for (int i = 0; i < npar; ++i)
        std::printf(" %.3f", std::sqrt(observedCovariance(i, i) / gnCovariance(i, i)));
      std::printf("\n");
    }
    positive = gn::covarianceFromHessian(hessian, covariance);
    if (!positive) hf_errlog(2026102255, "W: GaussNewton: Hessian not positive definite, pseudo-inverse used");
    _parameterCovariance.assign(npar * npar, 0.);
    for (int i = 0; i < npar; ++i)
      for (int j = 0; j < npar; ++j) _parameterCovariance[i * npar + j] = covariance(i, j);
  }
  const Eigen::VectorXd errors = covariance.diagonal().cwiseSqrt();
  _parameterUncertainties.assign(errors.data(), errors.data() + npar);
  for (int i = 0; i < npar; ++i)
    if (result.atBound[i]) hf_errlog(2026102256, "W: GaussNewton: parameter at a bound: " + _allParameterNames[i]);

  if (!result.converged) _convergenceStatus = ConvergenceStatus::NO_CONVERGENCE;
  else if (!positive) _convergenceStatus = ConvergenceStatus::FORCED_POSITIVE;
  else _convergenceStatus = ConvergenceStatus::SUCCESS;

  // Outputs, then the final fcn call (flag 3) at the minimum
  const gn::ProfiledChi2Diagnostics& diag = problem.diagnostics();
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  // the file summary must be reproducible (tools/test.sh): no timings, no round-off checks
  char line[256];
  std::vector<std::string> summary;
  std::snprintf(line, sizeof(line), "GaussNewton: %s, chi2 = %.6f, %d iterations", result.message.c_str(),
                result.chi2, result.iterations);
  summary.push_back(line);
  std::snprintf(line, sizeof(line), "fcn calls %d (Hessian %d), theory evaluations %d", evaluator.fcnCalls(),
                hessianCalls, jacobian.evaluations());
  summary.push_back(line);
  for (const std::string& s : summary) std::printf("%s\n", s.c_str());
  std::printf("wall time %.1f s (fcn %.1f, Jacobian %.1f, algebra %.1f)\n", wall, diag.fcnTime,
              diag.jacobianTime, diag.algebraTime);
  std::printf("checks: max |module chi2 - fcn| = %.2e, max |H_bb^-1 g_b| = %.2e\n", diag.maxChi2Mismatch,
              diag.maxShiftChange);

  gn::writeParameters(_allParameterNames, result.x, errors);
  if (settings.covariance) gn::writeCovariance(_allParameterNames, covariance, summary);
  evaluator.chi2(result.x, 3);
  gn::printParameters(_allParameterNames, result.x, errors);
}

void GaussNewtonMinimizer::errorAnalysis() {
  const YAML::Node node = XFITTER_PARS::rootNode["GaussNewton"]["doErrors"];
  if (!node) return;
  const std::string bandType = node.as<std::string>();
  if (bandType == "Hesse") {
    if (!hasParameterCovariance()) {
      hf_errlog(2026102257, "S: GaussNewton doErrors=Hesse requires covariance=1");
      return;
    }
    hf_errlog(12020506, "I: Calculation of symmetric error bands required");
    errbandssym_();
  } else if (bandType != "None") {
    hf_errlog(2026102258, "W: GaussNewton: unknown doErrors type " + bandType);
  }
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
