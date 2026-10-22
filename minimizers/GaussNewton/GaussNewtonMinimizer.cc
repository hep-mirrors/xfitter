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
#include "bridge/InputPdfs.h"
#include "bridge/TheoryEvaluator.h"
#include "core/IterativeHessian.h"
#include "core/LevenbergMarquardt.h"
#include "core/ObservedHessian.h"
#include "jacobian/FactorisedJacobian.h"
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
  int gn_ntot_();
  void gn_pumplin_member_(const int& j, const int& idir, const int& mpar, double* theoVars);
  void gn_pumplin_finish_(const int& mpar, double* theoVars);
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
  bool factorised = false;      ///< J_t = K F, Eq. (14)
  gn::FactorisedOptions factorisedOptions;
  bool covariance = true;
  bool observedHessian = true;  ///< covariance from second differences of the chi2 (as HESSE)
  bool centralHessian = true;   ///< central off-diagonal second differences
  double hessianDelta = 1e-3;   ///< chi2 change targeted by the second-difference steps
  bool checkDerivatives = false;
  bool pumplin = false;  ///< doErrors: Pumplin, iterative Hessian and eigenvector sets, Eq. (18)
  gn::IterativeHessianOptions iterative;
};

/// The GaussNewton block of parameters.yaml (README.md, "Steering")
Settings readSettings() {
  Settings s;
  const YAML::Node node = XFITTER_PARS::rootNode["GaussNewton"];
  if (!node.IsMap()) return s;
  const std::set<std::string> known = {"maxIterations", "tolerance", "lambda0", "jacobian", "fdScheme",
                                       "fdStep", "threads", "spanTolerance", "covariance", "hessian",
                                       "hessianScheme", "hessianDelta", "doErrors", "checkDerivatives",
                                       "pumplinIterations", "pumplinDelta", "pumplinTolerance"};
  for (const auto& item : node)
    if (!known.count(item.first.as<std::string>()))
      hf_errlog(2026102250, "W: GaussNewton: unknown steering key " + item.first.as<std::string>());

  s.lm.maxIterations = node["maxIterations"].as<int>(s.lm.maxIterations);
  s.lm.tolerance = node["tolerance"].as<double>(s.lm.tolerance);
  s.lm.lambda0 = node["lambda0"].as<double>(s.lm.lambda0);
  const std::string jacobian = node["jacobian"].as<std::string>("FiniteDifference");
  if (jacobian != "FiniteDifference" && jacobian != "Factorised")
    hf_errlog(2026102251, "F: GaussNewton: jacobian must be FiniteDifference or Factorised");
  s.factorised = jacobian == "Factorised";
  s.factorisedOptions.spanTolerance = node["spanTolerance"].as<double>(s.factorisedOptions.spanTolerance);
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
  const std::string hessianScheme = node["hessianScheme"].as<std::string>("central");
  if (hessianScheme != "forward" && hessianScheme != "central")
    hf_errlog(2026102260, "F: GaussNewton: hessianScheme must be forward or central");
  s.centralHessian = hessianScheme == "central";
  s.hessianDelta = node["hessianDelta"].as<double>(s.hessianDelta);
  s.checkDerivatives = node["checkDerivatives"].as<int>(0) != 0;
  s.pumplin = node["doErrors"].as<std::string>("None") == "Pumplin";
  s.iterative.maxIterations = node["pumplinIterations"].as<int>(s.iterative.maxIterations);
  s.iterative.delta = node["pumplinDelta"].as<double>(s.iterative.delta);
  s.iterative.tolerance = node["pumplinTolerance"].as<double>(s.iterative.tolerance);
  s.iterative.central = s.centralHessian;
  if (s.pumplin && !s.covariance) hf_errlog(2026102280, "F: GaussNewton: doErrors: Pumplin needs covariance: 1");
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

  gn::FiniteDifferenceJacobian fdJacobian(evaluator, steps, lower, upper, settings.fd);
  gn::InputPdfs inputPdfs(evaluator);
  gn::FactorisedJacobian factorisedJacobian(fdJacobian, inputPdfs, settings.factorisedOptions);
  gn::TheoryJacobian& jacobian = settings.factorised ? static_cast<gn::TheoryJacobian&>(factorisedJacobian)
                                                     : static_cast<gn::TheoryJacobian&>(fdJacobian);
  gn::ProfiledChi2 problem(evaluator, bridge, jacobian, lower, upper);
  problem.setWorkers(settings.fd.workers);  // the observed Hessian, Eq. (10)

  // First call (fcn flag 1) and start-point checks
  problem.setFcnFlag(1);
  problem.chi2(p0);
  if (bridge.polarisationChi2() != 0)
    hf_errlog(2026102253, "F: GaussNewton: polarisation shifts are not supported");
  const gn::DataSnapshot first = bridge.snapshot();
  factorisedJacobian.setPointScale(first.sigma);
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
    if (settings.pumplin) {
      // replaces the observed Hessian: its first iteration is one, with steps of Delta chi2 = delta
      const gn::IterativeHessianResult iterated =
          gn::iterativeHessian(problem, result.x, result.chi2, hessian, settings.iterative);
      hessianCalls = iterated.chi2Calls;
      std::printf("GaussNewton: Pumplin iterative Hessian, %d iterations, %d chi2 calls, max |H_y/2 - 1| = %.3f%s\n",
                  iterated.iterations, iterated.chi2Calls, iterated.deviation,
                  iterated.converged ? "" : " (not converged)");
      if (!iterated.converged)
        hf_errlog(2026102281, "W: GaussNewton: Pumplin iterative Hessian not converged, see pumplinIterations");
      if (iterated.droppedPairs > 0)
        hf_errlog(2026102285, "W: GaussNewton: Pumplin: some correlations not measured (invalid or far points)");
      if (iterated.lowestChi2 < result.chi2 - 1e-2)
        hf_errlog(2026102282, "W: GaussNewton: a lower chi2 was found around the minimum: not converged");
      hessian = iterated.hessian;
      _pumplinDisplacements = iterated.displacements;
      _pumplinCenter = result.x;
    } else if (settings.observedHessian) {
      const gn::ObservedHessianResult observed = gn::observedHessian(
          problem, result.x, result.chi2, hessian, settings.hessianDelta, settings.centralHessian);
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
  if (settings.factorised) {
    const gn::FactorisedStatistics& stats = factorisedJacobian.statistics();
    std::printf("Factorised Jacobian: %d columns from K F, %d by finite differences, %d explored directions, "
                "%d resets; exact Jacobian from iteration %d\n", stats.responseColumns, stats.fdColumns,
                factorisedJacobian.dimension(), stats.resets, result.exactModelAt);
  }
  std::printf("wall time %.1f s (fcn %.1f, Jacobian %.1f, algebra %.1f)\n", wall, diag.fcnTime,
              diag.jacobianTime, diag.algebraTime);
  std::printf("checks: max |module chi2 - fcn| = %.2e, max |H_bb^-1 g_b| = %.2e\n", diag.maxChi2Mismatch,
              diag.maxShiftChange);
  if (bridge.config().extraSystRescale)
    std::printf("checks: ExtraSystRescale, max |sigma rebuilt - sigma| / sigma = %.2e\n", diag.maxRescaleMismatch);

  gn::writeParameters(_allParameterNames, result.x, errors, steps);
  if (settings.covariance) gn::writeCovariance(_allParameterNames, covariance, summary);
  evaluator.chi2(result.x, 3);
  gn::printParameters(_allParameterNames, result.x, errors, steps);
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
  } else if (bandType == "Pumplin") {
    pumplinBands();
  } else if (bandType != "None") {
    hf_errlog(2026102258, "W: GaussNewton: unknown doErrors type " + bandType);
  }
}

void GaussNewtonMinimizer::pumplinBands() {
  const int npar = static_cast<int>(_pumplinCenter.size());
  if (npar == 0 || _pumplinDisplacements.cols() != npar) {
    hf_errlog(2026102283, "S: GaussNewton doErrors=Pumplin: no eigenvectors (covariance=0?)");
    return;
  }
  if (npar > 40) {
    hf_errlog(2026102284, "S: GaussNewton doErrors=Pumplin: at most 40 parameters (file names)");
    return;
  }
  hf_errlog(12020506, "I: Calculation of error bands required");
  double** pars = getPars();
  std::vector<double> theoVars(static_cast<size_t>(gn_ntot_()) * 2 * npar, 0.);
  for (int j = 0; j < npar; ++j)
    for (const int idir : {-1, 1}) {
      for (int i = 0; i < npar; ++i) *pars[i] = _pumplinCenter(i) + idir * _pumplinDisplacements(i, j);
      gn_pumplin_member_(j + 1, idir, npar, theoVars.data());
    }
  for (int i = 0; i < npar; ++i) *pars[i] = _pumplinCenter(i);
  delete[] pars;
  gn_pumplin_finish_(npar, theoVars.data());
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
