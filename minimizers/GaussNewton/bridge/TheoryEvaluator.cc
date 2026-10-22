#include "bridge/TheoryEvaluator.h"
#include "bridge/FortranBridge.h"
#include "BaseMinimizer.h"

#include <algorithm>
#include <cmath>

extern "C" {
  void fcn_(const int& npar, const double& dummy, double& chi2out, const double* pars,
            const int& iflag, const double& dummy2);
  void update_theory_iteration_();
}

namespace xfitter {
namespace gn {

namespace {
const size_t kFcnParameters = 200;  // fcn copies MNE = 200 values (endmini.inc)
}

TheoryEvaluator::TheoryEvaluator(const BaseMinimizer& minimizer, const FortranBridge& bridge)
    : _bridge(bridge), _priors(*minimizer.getPriors()) {
  double** pars = minimizer.getPars();
  _parameters.assign(pars, pars + minimizer.getNpars());
  delete[] pars;
}

void TheoryEvaluator::setParameters(const Eigen::VectorXd& p) {
  for (size_t i = 0; i < _parameters.size(); ++i) *_parameters[i] = p(i);
}

double TheoryEvaluator::chi2(const Eigen::VectorXd& p, int iflag) {
  setParameters(p);
  // fcn reads the parameters from gParameters; the array is only stored for MINUIT
  std::vector<double> pars(std::max(kFcnParameters, static_cast<size_t>(p.size())), 0.);
  std::copy(p.data(), p.data() + p.size(), pars.begin());
  fcn_(static_cast<int>(p.size()), 0., _lastFcn, pars.data(), iflag, 0.);
  ++_fcnCalls;
  return _lastFcn + priorChi2(p);
}

Eigen::VectorXd TheoryEvaluator::theory(const Eigen::VectorXd& p) {
  setParameters(p);
  update_theory_iteration_();
  ++_theoryCalls;
  return _bridge.theory();
}

double TheoryEvaluator::priorChi2(const Eigen::VectorXd& p) const {
  double chi2 = 0;
  for (size_t i = 0; i < _priors.size(); ++i)
    if (_priors[i]) chi2 += std::pow((p(i) - _priors[i][0]) / _priors[i][1], 2);
  return chi2;
}

void TheoryEvaluator::addPriors(const Eigen::VectorXd& p, QuadraticModel& model) const {
  for (size_t i = 0; i < _priors.size(); ++i) {
    if (!_priors[i]) continue;
    const double s2 = _priors[i][1] * _priors[i][1];
    model.gradient(i) += 2 * (p(i) - _priors[i][0]) / s2;
    model.hessian(i, i) += 2 / s2;
  }
}

}
}
