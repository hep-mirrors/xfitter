#include "jacobian/FiniteDifferenceJacobian.h"
#include "bridge/TheoryEvaluator.h"
#include "util/ForkPool.h"

#include <algorithm>
#include <cmath>

namespace xfitter {
namespace gn {

FiniteDifferenceJacobian::FiniteDifferenceJacobian(TheoryEvaluator& evaluator,
                                                   const Eigen::VectorXd& yamlSteps,
                                                   const Eigen::VectorXd& lower,
                                                   const Eigen::VectorXd& upper,
                                                   const FDOptions& options)
    : _evaluator(evaluator), _yamlSteps(yamlSteps), _lower(lower), _upper(upper), _options(options) {}

double FiniteDifferenceJacobian::step(const Eigen::VectorXd& p, int j) const {
  return _options.relativeStep * std::max(std::abs(p(j)), _yamlSteps(j));
}

bool FiniteDifferenceJacobian::isCentral(const Eigen::VectorXd& p, int j) const {
  const double h = step(p, j);
  return _options.central && p(j) + h <= _upper(j) && p(j) - h >= _lower(j);
}

Eigen::VectorXd FiniteDifferenceJacobian::column(const Eigen::VectorXd& p, const Eigen::VectorXd& t0, int j) {
  const double h = step(p, j);
  Eigen::VectorXd q = p;
  if (isCentral(p, j)) {
    q(j) = p(j) + h;
    const Eigen::VectorXd up = _evaluator.theory(q);
    q(j) = p(j) - h;
    const Eigen::VectorXd down = _evaluator.theory(q);
    return (up - down) / (2 * h);
  }
  // one-sided, away from the bound; shortened if neither side has room for h
  const double up = _upper(j) - p(j), down = p(j) - _lower(j);
  double signedStep = up >= h ? h : (down >= h ? -h : (up >= down ? up : -down));
  if (signedStep == 0) return Eigen::VectorXd::Zero(t0.size());
  q(j) = p(j) + signedStep;
  return (_evaluator.theory(q) - t0) / signedStep;
}

Eigen::MatrixXd FiniteDifferenceJacobian::compute(const Eigen::VectorXd& p, const Eigen::VectorXd& t0) {
  const int n = static_cast<int>(t0.size()), P = static_cast<int>(p.size());
  const std::vector<double> columns = ForkPool(_options.workers).run(
      P, n, [&](int j, double* out) { Eigen::Map<Eigen::VectorXd>(out, n) = column(p, t0, j); });
  for (int j = 0; j < P; ++j) _evaluations += isCentral(p, j) ? 2 : 1;
  return Eigen::Map<const Eigen::MatrixXd>(columns.data(), n, P);
}

}
}
