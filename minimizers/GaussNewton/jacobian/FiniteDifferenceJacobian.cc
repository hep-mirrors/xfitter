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

FiniteDifferenceJacobian::Displacement FiniteDifferenceJacobian::displacement(const Eigen::VectorXd& p, int j,
                                                                              bool central) const {
  const double h = step(p, j);
  const double up = _upper(j) - p(j), down = p(j) - _lower(j);
  if (central && up >= h && down >= h) return {h, -h};
  // one-sided, away from the bound; shortened if neither side has room for h
  return {up >= h ? h : (down >= h ? -h : (up >= down ? up : -down)), 0.};
}

Eigen::VectorXd FiniteDifferenceJacobian::column(const Eigen::VectorXd& p, const Eigen::VectorXd& t0, int j) {
  const Displacement d = displacement(p, j);
  if (d.up == d.down) return Eigen::VectorXd::Zero(t0.size());
  Eigen::VectorXd q = p;
  q(j) = p(j) + d.up;
  const Eigen::VectorXd up = _evaluator.theory(q);
  if (d.down == 0) return (up - t0) / d.up;
  q(j) = p(j) + d.down;
  return (up - _evaluator.theory(q)) / (d.up - d.down);
}

Eigen::MatrixXd FiniteDifferenceJacobian::columns(const Eigen::VectorXd& p, const Eigen::VectorXd& t0,
                                                  const std::vector<int>& which) {
  const int n = static_cast<int>(t0.size()), P = static_cast<int>(which.size());
  const std::vector<double> values = ForkPool(_options.workers).run(
      P, n, [&](int k, double* out) { Eigen::Map<Eigen::VectorXd>(out, n) = column(p, t0, which[k]); });
  for (const int j : which) {
    const Displacement d = displacement(p, j);
    _evaluations += d.up == d.down ? 0 : (d.down == 0 ? 1 : 2);
  }
  return Eigen::Map<const Eigen::MatrixXd>(values.data(), n, P);
}

Eigen::MatrixXd FiniteDifferenceJacobian::compute(const Eigen::VectorXd& p, const Eigen::VectorXd& t0) {
  std::vector<int> all(p.size());
  for (int j = 0; j < p.size(); ++j) all[j] = j;
  return columns(p, t0, all);
}

}
}
