#include "core/IterativeHessian.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace xfitter {
namespace gn {

namespace {

/// The problem in the coordinates y of x = x0 + U y, without bounds; records the lowest chi2
class Transformed : public Problem {
public:
  Transformed(Problem& problem, const Eigen::VectorXd& x0, const Eigen::MatrixXd& U, double& lowest)
      : _problem(problem), _x0(x0), _U(U), _lowest(lowest) {}
  int size() const override { return static_cast<int>(_U.cols()); }
  double chi2(const Eigen::VectorXd& y) override { return record(_problem.chi2(_x0 + _U * y)); }
  std::vector<double> chi2Batch(const std::vector<Eigen::VectorXd>& ys) override {
    std::vector<Eigen::VectorXd> xs;
    for (const Eigen::VectorXd& y : ys) xs.push_back(_x0 + _U * y);
    std::vector<double> values = _problem.chi2Batch(xs);
    for (const double v : values) record(v);
    return values;
  }
  QuadraticModel model(const Eigen::VectorXd&) override { return {}; }

private:
  double record(double v) {
    if (v < _lowest) _lowest = v;
    return v;
  }
  Problem& _problem;
  Eigen::VectorXd _x0;
  Eigen::MatrixXd _U;
  double& _lowest;
};

/// Step along one eigenvector (y units) and the chi2 at +-step
struct Direction {
  double step = 1, fPlus = 0, fMinus = 0;
  bool measured = false;  ///< a finite, positive sagitta was found
};

Eigen::VectorXd axis(int n, int j, double s) {
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
  y(j) = s;
  return y;
}

/// Steps with sagitta (f(+s) + f(-s))/2 - f0 within [delta/2, 2 delta], as MINUIT ITERATE
/// (Eq. 18): up to maxCycles batches, each for the directions not yet within the range
void findSteps(Transformed& problem, double f0, const IterativeHessianOptions& options,
               std::vector<Direction>& dir, int& calls) {
  const int n = static_cast<int>(dir.size());
  std::vector<int> active(n);
  for (int j = 0; j < n; ++j) {
    active[j] = j;
    dir[j] = Direction();
    dir[j].step = std::sqrt(options.delta);
  }
  for (int cycle = 0; cycle < options.maxCycles && !active.empty(); ++cycle) {
    std::vector<Eigen::VectorXd> points;
    for (const int j : active) {
      points.push_back(axis(n, j, dir[j].step));
      points.push_back(axis(n, j, -dir[j].step));
    }
    const std::vector<double> f = problem.chi2Batch(points);
    calls += static_cast<int>(points.size());
    std::vector<int> again;
    for (size_t a = 0; a < active.size(); ++a) {
      Direction& d = dir[active[a]];
      const double sagitta = 0.5 * (f[2 * a] + f[2 * a + 1]) - f0;
      if (std::isfinite(sagitta) && sagitta > 0) {
        d.fPlus = f[2 * a];
        d.fMinus = f[2 * a + 1];
        d.measured = true;
        if (sagitta >= options.delta / 2 && sagitta <= 2 * options.delta) continue;
        if (cycle + 1 == options.maxCycles) continue;  // keep the last measured step
        d.step *= std::min(10., std::max(0.1, std::sqrt(options.delta / sagitta)));
      } else {
        // invalid point (e.g. an invalid PDF): shorter; flat or concave: longer
        d.step *= std::isfinite(sagitta) ? 3. : 0.3;
        d.measured = false;
      }
      again.push_back(active[a]);
    }
    active = again;
  }
}


/// One side of an eigenvector in the calibration: the multiple t of the quadratic step with
/// Delta chi2 = delta, bracketed between lo (below delta) and hi (above, or an invalid point)
struct Side {
  double t = 1, tLo = 0, dLo = 0, tHi = 0, dHi = 0;
  bool done = false, failed = false;
};

/// As the final step of MINUIT ITERATE (its "alternative procedure"): along +-U_j, the exact
/// multiple t_+-, to calibrationAccuracy, with Delta chi2 = delta at t sqrt(delta) U_j, all sides
/// in one batch per cycle; then 1/t^2 = 1/t_+^2 + 1/t_-^2 - 1/(t_+ t_-), the curvature of a
/// quadratic + cubic through both crossings, and U_j -> t U_j. One crossing: its t; none: t = 1.
void calibrate(Transformed& problem, double f0, const IterativeHessianOptions& options,
               IterativeHessianResult& result) {
  const int n = static_cast<int>(result.displacements.cols());
  const double delta = options.delta, far = 100 * delta;
  std::vector<Side> sides(2 * n);
  for (int cycle = 0; cycle < options.calibrationCycles; ++cycle) {
    std::vector<int> active;
    std::vector<Eigen::VectorXd> points;
    for (int s = 0; s < 2 * n; ++s) {
      if (sides[s].done || sides[s].failed) continue;
      active.push_back(s);
      points.push_back(axis(n, s / 2, (s % 2 ? -1 : 1) * sides[s].t * std::sqrt(delta)));
    }
    if (active.empty()) break;
    const std::vector<double> f = problem.chi2Batch(points);
    result.chi2Calls += static_cast<int>(points.size());
    for (size_t a = 0; a < active.size(); ++a) {
      Side& side = sides[active[a]];
      const double rise = f[a] - f0;
      if (std::isfinite(rise) && std::abs(rise - delta) < options.calibrationAccuracy * delta) {
        side.done = true;
        continue;
      }
      if (std::isfinite(rise) && rise <= 0) {
        side.failed = true;  // below the minimum, as MINUIT: no crossing on this side
        continue;
      }
      const bool valid = std::isfinite(rise) && rise < far;
      if (valid && rise < delta) {
        side.tLo = side.t;
        side.dLo = rise;
      } else {
        side.tHi = side.t;
        side.dHi = valid ? rise : 0;  // 0: invalid or far, only bounds the step
      }
      if (side.tLo > 0 && side.tHi > 0 && side.dHi > 0) {
        // power law through the bracket, kept inside it
        const double p = std::log(side.dHi / side.dLo) / std::log(side.tHi / side.tLo);
        side.t = side.tLo * std::pow(delta / side.dLo, 1 / p);
        if (!(side.t > side.tLo && side.t < side.tHi)) side.t = std::sqrt(side.tLo * side.tHi);
      } else if (side.tLo > 0 && side.tHi > 0) {
        side.t = std::sqrt(side.tLo * side.tHi);
      } else if (valid) {
        side.t *= std::min(2., std::max(0.5, std::sqrt(delta / rise)));  // quadratic, at most 2x
      } else {
        side.t *= 0.5;
      }
    }
  }
  result.calibration = Eigen::VectorXd::Ones(n);
  for (int j = 0; j < n; ++j) {
    const Side &plus = sides[2 * j], &minus = sides[2 * j + 1];
    double t = 1;
    if (plus.done && minus.done)
      t = 1 / std::sqrt(1 / (plus.t * plus.t) + 1 / (minus.t * minus.t) - 1 / (plus.t * minus.t));
    else if (plus.done)
      t = plus.t;
    else if (minus.done)
      t = minus.t;
    else
      ++result.uncalibrated;
    result.calibration(j) = t;
    result.displacements.col(j) *= t;
  }
}

}  // namespace

bool eigenDisplacements(const Eigen::MatrixXd& hessian, Eigen::MatrixXd& displacements) {
  const int n = static_cast<int>(hessian.rows());
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(0.5 * (hessian + hessian.transpose()));
  const Eigen::VectorXd& values = eigen.eigenvalues();  // ascending
  const double floor = 1e-12 * std::max(values.cwiseAbs().maxCoeff(), std::numeric_limits<double>::min());
  bool positive = true;
  displacements.resize(n, n);
  for (int k = 0; k < n; ++k) {
    const int i = n - 1 - k;  // largest curvature first
    double value = values(i);
    if (!(value > floor)) {
      value = floor;
      positive = false;
    }
    Eigen::VectorXd v = eigen.eigenvectors().col(i);
    Eigen::Index largest;
    v.cwiseAbs().maxCoeff(&largest);
    if (v(largest) < 0) v = -v;
    displacements.col(k) = v * std::sqrt(2 / value);
  }
  return positive;
}

IterativeHessianResult iterativeHessian(Problem& problem, const Eigen::VectorXd& x, double chi2x,
                                        const Eigen::MatrixXd& startHessian,
                                        const IterativeHessianOptions& options) {
  const int n = static_cast<int>(x.size());
  IterativeHessianResult result;
  result.hessian = startHessian;
  result.lowestChi2 = chi2x;
  Eigen::MatrixXd best = startHessian;
  double bestDeviation = std::numeric_limits<double>::infinity();
  const double far = chi2x + 100 * options.delta;  // off-diagonal points above this are dropped

  for (int iteration = 0; iteration < options.maxIterations; ++iteration) {
    Eigen::MatrixXd U;
    eigenDisplacements(result.hessian, U);
    // chi2 = chi2x + |y|^2 if the current Hessian is right
    Transformed transformed(problem, x, U, result.lowestChi2);
    std::vector<Direction> dir(n);
    findSteps(transformed, chi2x, options, dir, result.chi2Calls);

    Eigen::MatrixXd Hy = Eigen::MatrixXd::Zero(n, n);
    for (int j = 0; j < n; ++j) {
      const Direction& d = dir[j];
      // a direction without a usable sagitta keeps the current curvature
      Hy(j, j) = d.measured ? (d.fPlus + d.fMinus - 2 * chi2x) / (d.step * d.step) : 2.;
    }
    // off-diagonal: as MINUIT ITERATE (corners), f at (+-a, +-b) with a, b the diagonal steps / sqrt(2),
    // which cancels the y_j^2 y_k^2 term and halves the y_j^3 y_k ones; otherwise f(++) [and f(--)]
    // with the diagonal steps
    std::vector<Eigen::VectorXd> points;
    std::vector<std::pair<int, int>> pairs;
    const double shrink = options.corners ? 1 / std::sqrt(2.) : 1.;
    for (int j = 0; j < n; ++j)
      for (int k = j + 1; k < n; ++k) {
        if (!dir[j].measured || !dir[k].measured) continue;
        pairs.push_back({j, k});
        const Eigen::VectorXd a = axis(n, j, shrink * dir[j].step), b = axis(n, k, shrink * dir[k].step);
        points.push_back(a + b);
        if (options.corners || options.central) points.push_back(-a - b);
        if (options.corners) {
          points.push_back(a - b);
          points.push_back(-a + b);
        }
      }
    const std::vector<double> f = transformed.chi2Batch(points);
    result.chi2Calls += static_cast<int>(points.size());
    const size_t perPair = options.corners ? 4 : options.central ? 2 : 1;
    for (size_t a = 0; a < pairs.size(); ++a) {
      const int j = pairs[a].first, k = pairs[a].second;
      const Direction &dj = dir[j], &dk = dir[k];
      const double* fa = f.data() + perPair * a;
      if (!std::all_of(fa, fa + perPair, [far](double v) { return v < far; })) {
        ++result.droppedPairs;  // an invalid or far point: no correlation measured
        continue;
      }
      if (options.corners)
        Hy(j, k) = (fa[0] + fa[1] - fa[2] - fa[3]) / (4 * shrink * shrink * dj.step * dk.step);
      else if (options.central)
        Hy(j, k) = (fa[0] + fa[1] - dj.fPlus - dj.fMinus - dk.fPlus - dk.fMinus + 2 * chi2x) / (2 * dj.step * dk.step);
      else
        Hy(j, k) = (fa[0] - dj.fPlus - dk.fPlus + chi2x) / (dj.step * dk.step);
      Hy(k, j) = Hy(j, k);
    }
    ++result.iterations;
    const double deviation = (0.5 * Hy - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff();

    // back to x: H = U^-T H_y U^-1, kept positive definite (errors grow at most 10x per iteration)
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(0.5 * (Hy + Hy.transpose()));
    const Eigen::VectorXd values = eigen.eigenvalues().cwiseMax(0.02);
    const Eigen::MatrixXd HyPositive = eigen.eigenvectors() * values.asDiagonal() * eigen.eigenvectors().transpose();
    const Eigen::MatrixXd Uinv = U.inverse();
    const Eigen::MatrixXd H = Uinv.transpose() * HyPositive * Uinv;
    result.hessian = 0.5 * (H + H.transpose());
    if (deviation < bestDeviation) {
      bestDeviation = deviation;
      best = result.hessian;
    }
    if (deviation < options.tolerance) {
      result.converged = true;
      break;
    }
  }
  if (!result.converged) result.hessian = best;
  result.deviation = bestDeviation;
  eigenDisplacements(result.hessian, result.displacements);
  if (options.calibrate) {
    Transformed transformed(problem, x, result.displacements, result.lowestChi2);
    calibrate(transformed, chi2x, options, result);
  }
  return result;
}

}
}
