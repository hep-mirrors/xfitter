// Standalone tests of minimizers/GaussNewton/core (no xFitter needed).
// Build with -DXFITTER_GN_TESTS=ON and run gn_core_tests; exit code 0 means all passed.

#include "core/LevenbergMarquardt.h"
#include "core/ProfiledNuisances.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using namespace xfitter::gn;

namespace {

int failures = 0;

std::string sci(double x) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.1e", x);
  return buffer;
}

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


/// Toy version of the xFitter chi2, Eq. (1): t(p) = T p + c (p_0 p_1)^2 (mildly nonlinear),
/// theory-dependent sigma and Gamma, log term, nuisance shifts profiled exactly.
struct ToyChi2 {
  int n = 30, m = 5, P = 4;
  bool linear = false;  // no t-dependence of sigma and Gamma, no log term, linear t(p)
  Eigen::MatrixXd T, beta;
  Eigen::VectorXd d, e_dt, e_tt, e_const, power, prior, alpha;
  std::vector<bool> additive, logTerm;

  explicit ToyChi2(bool lin) : linear(lin) {
    std::srand(7);
    T = (Eigen::MatrixXd::Random(n, P).array() * 0.3 + 1).matrix();
    beta = Eigen::MatrixXd::Random(m, n) * 0.03;
    d = Eigen::VectorXd::Constant(n, 2.0) + 0.2 * Eigen::VectorXd::Random(n);
    e_dt = Eigen::VectorXd::Constant(n, 0.02 * 0.02);
    e_tt = Eigen::VectorXd::Constant(n, 0.01 * 0.01);
    for (int i = 0; i < n; i += 2) e_dt(i) = 0;  // these points have d2L/dt2 > 0
    e_const = Eigen::VectorXd::Constant(n, 0.015 * 0.015);
    power = Eigen::VectorXd(m);
    power << 1, 0.5, 0, 1, 1;
    prior = Eigen::VectorXd::Constant(m, 1.0);
    prior(2) = 0.5;
    alpha = Eigen::VectorXd::Constant(n, 0.05);
    additive.assign(n, false);
    additive[3] = true;
    logTerm.assign(n, !linear);
    if (linear) {
      e_dt.setZero();
      e_tt.setZero();
      power.setZero();
    }
  }

  Eigen::VectorXd theory(const Eigen::VectorXd& p) const {
    Eigen::VectorXd t = T * p;
    if (!linear) t.array() += 0.05 * std::pow(p(0) * p(1), 2);
    return t;
  }

  /// Snapshot at theory t and shifts b, computed the way xFitter does
  DataSnapshot snapshot(const Eigen::VectorXd& t, const Eigen::VectorXd& b) const {
    DataSnapshot s;
    s.data = d;
    s.theory = t;
    s.sigma.resize(n);
    s.sigma2DT = e_dt;
    s.sigma2TT = e_tt;
    s.additive = additive;
    s.logTerm = logTerm;
    s.alpha = alpha;
    s.gamma.resize(m, n);
    s.gammaPower = power;
    s.shifts = b;
    s.priorScale = prior;
    for (int i = 0; i < n; ++i) {
      const double tt = t(i) > 0 ? t(i) : d(i);
      s.sigma(i) = std::sqrt(e_dt(i) * std::abs(d(i) * tt) + e_tt(i) * tt * tt + e_const(i) * d(i) * d(i));
      for (int k = 0; k < m; ++k) {
        const double scale = additive[i] ? d(i) : std::pow(t(i), power(k)) * std::pow(d(i), 1 - power(k));
        s.gamma(k, i) = beta(k, i) * scale;
      }
    }
    return s;
  }

  /// Exact profile of b for fixed t: r is linear in b
  Eigen::VectorXd profile(const Eigen::VectorXd& t) const {
    const DataSnapshot s = snapshot(t, Eigen::VectorXd::Zero(m));
    const PointModel pm = computePointModel(s);
    Eigen::MatrixXd A = pm.dResidualdB * pm.dResidualdB.transpose();
    A.diagonal() += prior;
    return A.ldlt().solve(-pm.dResidualdB * pm.residual);
  }

  double chi2(const Eigen::VectorXd& p) const {
    const Eigen::VectorXd t = theory(p);
    return computePointModel(snapshot(t, profile(t))).chi2();
  }

  Eigen::MatrixXd jacobian(const Eigen::VectorXd& p) const {
    Eigen::MatrixXd J(n, P);
    for (int j = 0; j < P; ++j) {
      Eigen::VectorXd up = p, dn = p;
      up(j) += 1e-6;
      dn(j) -= 1e-6;
      J.col(j) = (theory(up) - theory(dn)) / 2e-6;
    }
    return J;
  }
};

Eigen::VectorXd toyPoint() {
  Eigen::VectorXd p(4);
  p << 0.7, 0.9, 0.4, 0.5;
  return p;
}

void testPointDerivatives() {
  const ToyChi2 toy(false);
  const Eigen::VectorXd t = toy.theory(toyPoint());
  const Eigen::VectorXd b = toy.profile(t);
  const PointModel pm = computePointModel(toy.snapshot(t, b));
  double worst = 0, worstLog = 0, worstCurv = 0;
  int positiveCurvature = 0;
  for (int i = 0; i < toy.n; ++i) {
    const double h = 1e-5 * t(i);
    Eigen::VectorXd up = t, dn = t;
    up(i) += h;
    dn(i) -= h;
    const PointModel pu = computePointModel(toy.snapshot(up, b));
    const PointModel pd = computePointModel(toy.snapshot(dn, b));
    const double a = (pu.residual(i) - pd.residual(i)) / (2 * h);
    const double l1 = (pu.logTerm(i) - pd.logTerm(i)) / (2 * h);
    const double l2 = (pu.logTerm(i) - 2 * pm.logTerm(i) + pd.logTerm(i)) / (h * h);
    worst = std::max(worst, std::abs(a - pm.dResidualdT(i)) / std::abs(a));
    worstLog = std::max(worstLog, std::abs(l1 - pm.dLogTermdT(i)) / std::abs(l1));
    worstCurv = std::max(worstCurv, std::abs(std::max(l2, 0.) - pm.d2LogTermdT2(i)) / (std::abs(l2) + 1e-3));
    if (l2 > 0) ++positiveCurvature;
  }
  check(positiveCurvature > 0 && positiveCurvature < toy.n, "log-term curvature has both signs");
  check(worst < 1e-6, "dr/dt, Eq. (5), matches finite differences: " + sci(worst));
  check(worstLog < 1e-6, "dL/dt, Eq. (6), matches finite differences: " + sci(worstLog));
  check(worstCurv < 1e-3, "max(d2L/dt2, 0), Eq. (6), matches finite differences: " + sci(worstCurv));
}

void testProfiledGradient() {
  const ToyChi2 toy(false);
  const Eigen::VectorXd p = toyPoint();
  const Eigen::VectorXd t = toy.theory(p);
  const DataSnapshot s = toy.snapshot(t, toy.profile(t));
  const PointModel pm = computePointModel(s);
  const ProfiledModel prof = profileNuisances(pm, toy.jacobian(p), s.priorScale);
  check(std::abs(pm.chi2() - toy.chi2(p)) < 1e-12 * toy.chi2(p), "module chi2 equals the toy chi2");
  check(prof.shiftGradient.norm() < 1e-9, "g_b = 0 at the exact profile, Eq. (9)");
  check(prof.shiftCorrection.norm() < 1e-9, "model shift correction = 0 at the exact profile");
  double worst = 0;
  for (int j = 0; j < toy.P; ++j) {
    Eigen::VectorXd up = p, dn = p;
    up(j) += 1e-5;
    dn(j) -= 1e-5;
    const double fd = (toy.chi2(up) - toy.chi2(dn)) / 2e-5;
    worst = std::max(worst, std::abs(fd - prof.model.gradient(j)) / std::abs(fd));
  }
  check(worst < 1e-5, "profiled gradient, Eq. (9), matches FD of the profiled chi2: " + sci(worst));
}

void testSchurVersusJoint() {
  // Linear case: the profiled Hessian is exact, and the profiled Newton step equals the
  // p-part of the joint (p, b) Newton step.
  const ToyChi2 toy(true);
  const Eigen::VectorXd p = toyPoint();
  const Eigen::VectorXd t = toy.theory(p);
  const Eigen::VectorXd b = Eigen::VectorXd::Constant(toy.m, 0.3);  // not the profile
  const DataSnapshot s = toy.snapshot(t, b);
  const PointModel pm = computePointModel(s);
  const Eigen::MatrixXd J = toy.jacobian(p);
  const ProfiledModel prof = profileNuisances(pm, J, s.priorScale);

  // joint residual vector and Jacobian, Eq. (7)
  const int n = toy.n, m = toy.m, P = toy.P;
  Eigen::MatrixXd Jj = Eigen::MatrixXd::Zero(n + m, P + m);
  Jj.topLeftCorner(n, P) = pm.dResidualdT.asDiagonal() * J;
  Jj.topRightCorner(n, m) = pm.dResidualdB.transpose();
  Jj.bottomRightCorner(m, m) = s.priorScale.cwiseSqrt().asDiagonal();
  Eigen::VectorXd rj(n + m);
  rj << pm.residual, pm.prior;
  const Eigen::VectorXd joint = (Jj.transpose() * Jj).ldlt().solve(-Jj.transpose() * rj);
  const Eigen::VectorXd schur = prof.model.hessian.ldlt().solve(-prof.model.gradient);
  check((joint.head(P) - schur).norm() < 1e-9 * joint.head(P).norm(), "Schur step equals the joint step");
  check((b + joint.tail(m) - toy.profile(toy.theory(p + schur))).norm() < 1e-9,
        "joint step lands on the exact profile");

  // the profiled Hessian equals the FD Hessian of the profiled chi2
  double worst = 0;
  for (int j = 0; j < P; ++j)
    for (int k = 0; k < P; ++k) {
      const double h = 1e-3;
      Eigen::VectorXd pp = p, pm_ = p, mp = p, mm = p;
      pp(j) += h; pp(k) += h;
      pm_(j) += h; pm_(k) -= h;
      mp(j) -= h; mp(k) += h;
      mm(j) -= h; mm(k) -= h;
      const double fd = (toy.chi2(pp) - toy.chi2(pm_) - toy.chi2(mp) + toy.chi2(mm)) / (4 * h * h);
      worst = std::max(worst, std::abs(fd - prof.model.hessian(j, k)) / prof.model.hessian.norm());
    }
  check(worst < 1e-6, "profiled Hessian, Eq. (9), matches FD in the linear case: " + sci(worst));
}

void testCovarianceFallback() {
  Eigen::MatrixXd H(2, 2);
  H << 2, 0, 0, 0;  // second direction flat
  Eigen::MatrixXd C;
  const bool pd = covarianceFromHessian(H, C);
  check(!pd && std::abs(C(0, 0) - 1) < 1e-12 && C(1, 1) == 0, "covariance pseudo-inverse fallback, Eq. (10)");
}

}  // namespace

int main() {
  testLinearLeastSquares();
  testRosenbrock();
  testBounds();
  testCorrelatedBounds();
  testCovarianceScale();
  testPointDerivatives();
  testProfiledGradient();
  testSchurVersusJoint();
  testCovarianceFallback();
  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
  return failures ? 1 : 0;
}
