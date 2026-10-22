// Standalone tests of minimizers/GaussNewton/core (no xFitter needed).
// Build with -DXFITTER_GN_TESTS=ON and run gn_core_tests; exit code 0 means all passed.

#include "core/LevenbergMarquardt.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using namespace xfitter::gn;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  std::printf("%s  %s\n", ok ? "PASS" : "FAIL", what.c_str());
  if (!ok) ++failures;
}

/// chi2 = |r(x)|^2 with Gauss-Newton model g = 2 J^T r, H = 2 J^T J
class LeastSquares : public Problem {
public:
  using Residuals = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;
  using Jacobian = std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>;
  LeastSquares(int n, Residuals r, Jacobian J) : _n(n), _r(r), _J(J) {}
  int size() const override { return _n; }
  double chi2(const Eigen::VectorXd& x) override { return _r(x).squaredNorm(); }
  QuadraticModel model(const Eigen::VectorXd& x) override {
    const Eigen::MatrixXd J = _J(x);
    return {2 * J.transpose() * _r(x), 2 * J.transpose() * J};
  }
  Eigen::VectorXd lower, upper;
  Eigen::VectorXd lowerBounds() const override { return lower.size() ? lower : Problem::lowerBounds(); }
  Eigen::VectorXd upperBounds() const override { return upper.size() ? upper : Problem::upperBounds(); }

private:
  int _n;
  Residuals _r;
  Jacobian _J;
};

// tolerance 1e-12 in chi2 units: x is then known to ~1e-6 / sqrt(curvature)
LMOptions tightOptions() {
  LMOptions options;
  options.tolerance = 1e-12;
  options.maxIterations = 200;
  return options;
}

void testLinearLeastSquares() {
  Eigen::MatrixXd A(5, 3);
  A << 1, 2, 0, 0, 1, 1, 3, 0, 1, 1, 1, 1, 2, -1, 0;
  Eigen::VectorXd y(5);
  y << 1, 2, 0, -1, 3;
  LeastSquares problem(3, [&](const Eigen::VectorXd& x) { return Eigen::VectorXd(A * x - y); },
                       [&](const Eigen::VectorXd&) { return A; });
  const LMResult result = LevenbergMarquardt(tightOptions()).minimise(problem, Eigen::VectorXd::Zero(3));
  const Eigen::VectorXd exact = (A.transpose() * A).ldlt().solve(A.transpose() * y);
  check(result.converged, "linear LSQ converges: " + result.message);
  check((result.x - exact).norm() < 1e-8, "linear LSQ matches the normal equations");
  check(result.iterations <= 3, "linear LSQ needs at most 3 iterations, got " +
                                    std::to_string(result.iterations));
}

void testRosenbrock() {
  LeastSquares problem(
      2,
      [](const Eigen::VectorXd& x) {
        Eigen::VectorXd r(2);
        r << 10 * (x(1) - x(0) * x(0)), 1 - x(0);
        return r;
      },
      [](const Eigen::VectorXd& x) {
        Eigen::MatrixXd J(2, 2);
        J << -20 * x(0), 10, -1, 0;
        return J;
      });
  Eigen::VectorXd x0(2);
  x0 << -1.2, 1;
  const LMResult result = LevenbergMarquardt(tightOptions()).minimise(problem, x0);
  check(result.converged, "Rosenbrock converges: " + result.message);
  check((result.x - Eigen::Vector2d(1, 1)).norm() < 1e-6, "Rosenbrock minimum at (1,1)");
  std::printf("      Rosenbrock: %d iterations, %d chi2 calls\n", result.iterations, result.chi2Calls);
}

void testBounds() {
  // chi2 = sum_i w_i (x_i - c_i)^2, separable, so the bounded minimum is the clamped c
  Eigen::VectorXd c(3), w(3);
  c << 2, -3, 0.5;
  w << 1, 4, 9;
  LeastSquares problem(3,
                       [&](const Eigen::VectorXd& x) {
                         return Eigen::VectorXd(w.cwiseSqrt().cwiseProduct(x - c));
                       },
                       [&](const Eigen::VectorXd&) { return Eigen::MatrixXd(w.cwiseSqrt().asDiagonal()); });
  problem.lower = Eigen::Vector3d(-1, -1, -1);
  problem.upper = Eigen::Vector3d(1, 1, 1);
  const LMResult result = LevenbergMarquardt(tightOptions()).minimise(problem, Eigen::VectorXd::Zero(3));
  check(result.converged, "bounded problem converges: " + result.message);
  check((result.x - Eigen::Vector3d(1, -1, 0.5)).norm() < 1e-6, "bounded minimum is the clamped one");
  check(result.atBound[0] && result.atBound[1] && !result.atBound[2], "bound flags");
}

void testCorrelatedBounds() {
  // Non-separable quadratic with one active bound: compare with the reduced problem
  Eigen::MatrixXd H(2, 2);
  H << 4, 1, 1, 2;
  Eigen::VectorXd c(2);
  c << 2, 1;  // unconstrained minimum
  const Eigen::LLT<Eigen::MatrixXd> llt(H);
  const Eigen::MatrixXd L = llt.matrixU();
  LeastSquares problem(2, [&](const Eigen::VectorXd& x) { return Eigen::VectorXd(L * (x - c)); },
                       [&](const Eigen::VectorXd&) { return L; });
  problem.lower = Eigen::Vector2d(-5, -5);
  problem.upper = Eigen::Vector2d(1, 5);
  const LMResult result = LevenbergMarquardt(tightOptions()).minimise(problem, Eigen::VectorXd::Zero(2));
  // x0 = 1 at the bound; minimise over x1: H11 (x1 - c1) + H10 (x0 - c0) = 0
  const double x1 = c(1) - H(1, 0) * (1 - c(0)) / H(1, 1);
  check(result.converged, "correlated bounded problem converges: " + result.message);
  check(std::abs(result.x(0) - 1) < 1e-10 && std::abs(result.x(1) - x1) < 1e-6,
        "correlated bounded minimum");
}

void testCovarianceScale() {
  // For chi2 = sum ((x - y_i)/s)^2 with n points, the error on x is s/sqrt(n): C = 2 H^-1
  const int n = 4;
  const double s = 0.5;
  LeastSquares problem(1,
                       [&](const Eigen::VectorXd& x) {
                         Eigen::VectorXd r(n);
                         for (int i = 0; i < n; ++i) r(i) = (x(0) - i) / s;
                         return r;
                       },
                       [&](const Eigen::VectorXd&) { return Eigen::MatrixXd::Constant(n, 1, 1 / s); });
  const LMResult result = LevenbergMarquardt(tightOptions()).minimise(problem, Eigen::VectorXd::Zero(1));
  const double error = std::sqrt(2 / result.model.hessian(0, 0));
  check(std::abs(error - s / std::sqrt(n)) < 1e-12, "C = 2 H^-1 gives the Delta chi2 = 1 error");
}

}  // namespace

int main() {
  testLinearLeastSquares();
  testRosenbrock();
  testBounds();
  testCorrelatedBounds();
  testCovarianceScale();
  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
  return failures ? 1 : 0;
}
