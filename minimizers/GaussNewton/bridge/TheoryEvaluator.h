#pragma once

// The two ways the module evaluates xFitter (README.md, Sec. 2.3):
// - chi2(p): the official chi2, fcn plus the Gaussian parameter priors, Eq. (1);
// - theory(p): theory only, through update_theory_iteration, for Jacobian columns.
// Both change the global xFitter state (parameters, PDFs, THEO).

#include "core/Problem.h"

#include <vector>

namespace xfitter {
class BaseMinimizer;

namespace gn {
class FortranBridge;

class TheoryEvaluator {
public:
  TheoryEvaluator(const BaseMinimizer& minimizer, const FortranBridge& bridge);

  /// fcn(iflag) + sum ((p - mu)/s)^2. iflag: 1 first call, 2 iteration, 3 final.
  double chi2(const Eigen::VectorXd& p, int iflag = 2);

  /// fcn value of the last chi2 call, without the parameter priors
  double lastFcn() const { return _lastFcn; }

  /// Theory of the fitted points at p, without computing the chi2
  Eigen::VectorXd theory(const Eigen::VectorXd& p);

  /// Gaussian parameter priors: value, and their gradient and Hessian added to a model
  double priorChi2(const Eigen::VectorXd& p) const;
  void addPriors(const Eigen::VectorXd& p, QuadraticModel& model) const;

  int fcnCalls() const { return _fcnCalls; }
  /// Adds fcn calls made by forked worker processes
  void countFcnCalls(int n) { _fcnCalls += n; }
  int theoryCalls() const { return _theoryCalls; }

  /// Writes p into XFITTER_PARS::gParameters; nothing is recomputed
  void setParameters(const Eigen::VectorXd& p);

private:
  const FortranBridge& _bridge;
  std::vector<double*> _parameters;     ///< pointers into XFITTER_PARS::gParameters
  std::vector<const double*> _priors;   ///< {mean, sigma} or nullptr, per parameter
  double _lastFcn = 0;
  int _fcnCalls = 0;
  int _theoryCalls = 0;
};

}
}
