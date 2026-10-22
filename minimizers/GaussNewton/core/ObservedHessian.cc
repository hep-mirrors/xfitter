#include "core/ObservedHessian.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace xfitter {
namespace gn {

namespace {

/// Diagonal scheme along e_j: central (x +- h), or one-sided away from a bound (x + s h, x + 2 s h)
struct Diagonal {
  double c = 0;  ///< curvature that sets the step
  double h = 0, side = 1;
  bool central = true;
  double fSide = 0, fOther = 0;  ///< chi2 at x + side h and at the other point of the scheme
  double measured = 0;
};

}  // namespace

ObservedHessianResult observedHessian(Problem& problem, const Eigen::VectorXd& x, double f0,
                                      const Eigen::MatrixXd& modelHessian, double targetDelta,
                                      bool centralOffDiagonal) {
  const int n = static_cast<int>(x.size());
  const Eigen::VectorXd lower = problem.lowerBounds(), upper = problem.upperBounds();
  ObservedHessianResult result;
  result.hessian = Eigen::MatrixXd::Zero(n, n);
  std::vector<Diagonal> diag(n);
  const double tiny = 1e-12 * std::max(1., modelHessian.diagonal().cwiseAbs().maxCoeff());
  auto shifted = [&](int j, double hj, int k, double hk) {
    Eigen::VectorXd y = x;
    y(j) += hj;
    if (k >= 0) y(k) += hk;
    return y;
  };

  // Diagonal: one batch of points per pass; a second pass for the steps far from the
  // measured curvature scale
  std::vector<int> active(n);
  for (int j = 0; j < n; ++j) {
    active[j] = j;
    diag[j].c = std::max(modelHessian(j, j), tiny);
  }
  for (int pass = 0; pass < 2 && !active.empty(); ++pass) {
    std::vector<Eigen::VectorXd> points;
    for (const int j : active) {
      Diagonal& d = diag[j];
      // chi2 changes by 1/2 c h^2 = targetDelta for curvature c
      d.h = std::sqrt(2 * targetDelta / d.c);
      // stay inside the bounds: the one-sided scheme needs 2h on one side
      const double room = std::max(upper(j) - x(j), x(j) - lower(j));
      if (room > 0) d.h = std::min(d.h, room / 2);
      const bool upOk = x(j) + 2 * d.h <= upper(j), downOk = x(j) - 2 * d.h >= lower(j);
      d.central = upOk && downOk;
      d.side = d.central || upOk ? 1 : -1;
      points.push_back(shifted(j, d.side * d.h, -1, 0));
      points.push_back(shifted(j, d.central ? -d.h : 2 * d.side * d.h, -1, 0));
    }
    const std::vector<double> f = problem.chi2Batch(points);
    result.chi2Calls += static_cast<int>(points.size());
    std::vector<int> redo;
    for (size_t a = 0; a < active.size(); ++a) {
      Diagonal& d = diag[active[a]];
      d.fSide = f[2 * a];
      d.fOther = f[2 * a + 1];
      // central: (f(h) + f(-h) - 2 f0) / h^2; one-sided: (f(2sh) - 2 f(sh) + f0) / h^2
      d.measured = d.central ? (d.fSide + d.fOther - 2 * f0) / (d.h * d.h)
                             : (d.fOther - 2 * d.fSide + f0) / (d.h * d.h);
      result.hessian(active[a], active[a]) = d.measured;
      if (d.measured > 0 && !(d.measured < 4 * d.c && d.measured > d.c / 4)) {
        d.c = d.measured;  // redo once with the measured curvature
        redo.push_back(active[a]);
      }
    }
    active = redo;
  }

  // Off-diagonal: all points in one batch
  std::vector<Eigen::VectorXd> points;
  for (int j = 0; j < n; ++j)
    for (int k = j + 1; k < n; ++k) {
      points.push_back(shifted(j, diag[j].side * diag[j].h, k, diag[k].side * diag[k].h));
      if (centralOffDiagonal && diag[j].central && diag[k].central)
        points.push_back(shifted(j, -diag[j].h, k, -diag[k].h));
    }
  const std::vector<double> f = problem.chi2Batch(points);
  result.chi2Calls += static_cast<int>(points.size());
  size_t next = 0;
  for (int j = 0; j < n; ++j)
    for (int k = j + 1; k < n; ++k) {
      const Diagonal &dj = diag[j], &dk = diag[k];
      const double fjk = f[next++];
      if (centralOffDiagonal && dj.central && dk.central) {
        // f(++) + f(--) - f(+j) - f(-j) - f(+k) - f(-k) + 2 f0 = 2 H_jk h_j h_k + O(h^4)
        const double fmm = f[next++];
        result.hessian(j, k) = result.hessian(k, j) =
            (fjk + fmm - dj.fSide - dj.fOther - dk.fSide - dk.fOther + 2 * f0) / (2 * dj.h * dk.h);
      } else {
        result.hessian(j, k) = result.hessian(k, j) =
            (fjk - dj.fSide - dk.fSide + f0) / (dj.side * dj.h * dk.side * dk.h);
      }
    }
  return result;
}

}
}
