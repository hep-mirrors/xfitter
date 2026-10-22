#include "core/ObservedHessian.h"

#include <algorithm>
#include <cmath>

namespace xfitter {
namespace gn {

namespace {

/// Diagonal second derivative along e_j with step h (central, or one-sided at a bound).
/// Returns the curvature; side and fSide are the point used for the off-diagonal terms.
double diagonal(Problem& problem, const Eigen::VectorXd& x, double f0, int j, double h,
                const Eigen::VectorXd& lower, const Eigen::VectorXd& upper, int& calls,
                double& side, double& fSide) {
  auto f = [&](double shift) {
    Eigen::VectorXd y = x;
    y(j) += shift;
    ++calls;
    return problem.chi2(y);
  };
  const bool upOk = x(j) + 2 * h <= upper(j), downOk = x(j) - 2 * h >= lower(j);
  if (upOk && downOk) {
    side = 1;
    fSide = f(h);
    return (fSide + f(-h) - 2 * f0) / (h * h);
  }
  side = upOk ? 1 : -1;  // one-sided, away from the bound
  fSide = f(side * h);
  return (f(2 * side * h) - 2 * fSide + f0) / (h * h);
}

}  // namespace

ObservedHessianResult observedHessian(Problem& problem, const Eigen::VectorXd& x, double f0,
                                      const Eigen::MatrixXd& modelHessian, double targetDelta) {
  const int n = static_cast<int>(x.size());
  const Eigen::VectorXd lower = problem.lowerBounds(), upper = problem.upperBounds();
  ObservedHessianResult result;
  result.hessian = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd h(n), side(n), fSide(n);
  const double tiny = 1e-12 * std::max(1., modelHessian.diagonal().cwiseAbs().maxCoeff());

  for (int j = 0; j < n; ++j) {
    // chi2 changes by 1/2 c h^2 = targetDelta for curvature c
    double c = std::max(modelHessian(j, j), tiny);
    for (int pass = 0; pass < 2; ++pass) {
      h(j) = std::sqrt(2 * targetDelta / c);
      const double measured = diagonal(problem, x, f0, j, h(j), lower, upper, result.chi2Calls, side(j), fSide(j));
      result.hessian(j, j) = measured;
      if (!(measured > 0) || (measured < 4 * c && measured > c / 4)) break;
      c = measured;  // step far from the measured curvature scale: redo once
    }
  }

  for (int j = 0; j < n; ++j)
    for (int k = j + 1; k < n; ++k) {
      Eigen::VectorXd y = x;
      y(j) += side(j) * h(j);
      y(k) += side(k) * h(k);
      ++result.chi2Calls;
      const double fjk = problem.chi2(y);
      result.hessian(j, k) = result.hessian(k, j) =
          (fjk - fSide(j) - fSide(k) + f0) / (side(j) * h(j) * side(k) * h(k));
    }
  return result;
}

}
}
