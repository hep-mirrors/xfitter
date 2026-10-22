#include "core/LevenbergMarquardt.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace xfitter {
namespace gn {

namespace {

const double kAcceptRatio = 1e-4;  // Eq. (11)
const double kMinLambda = 1e-12;
const double kDampingFloor = 1e-8;  // relative floor of diag|H| in the damping, Eq. (11)

Eigen::VectorXd clampToBox(const Eigen::VectorXd& x, const Eigen::VectorXd& lower,
                           const Eigen::VectorXd& upper) {
  return x.cwiseMax(lower).cwiseMin(upper);
}

std::vector<int> indices(const std::vector<bool>& free) {
  std::vector<int> idx;
  for (size_t i = 0; i < free.size(); ++i)
    if (free[i]) idx.push_back(static_cast<int>(i));
  return idx;
}

/// Restriction of the model to the free parameters
void restrict(const QuadraticModel& model, const std::vector<int>& idx, Eigen::VectorXd& g,
              Eigen::MatrixXd& H) {
  const int n = static_cast<int>(idx.size());
  g.resize(n);
  H.resize(n, n);
  for (int i = 0; i < n; ++i) {
    g(i) = model.gradient(idx[i]);
    for (int j = 0; j < n; ++j) H(i, j) = model.hessian(idx[i], idx[j]);
  }
}

}  // namespace

std::vector<bool> freeParameters(const Eigen::VectorXd& x, const Eigen::VectorXd& gradient,
                                 const Eigen::VectorXd& lower, const Eigen::VectorXd& upper) {
  std::vector<bool> free(x.size(), true);
  for (int i = 0; i < x.size(); ++i) {
    // descent direction is -gradient
    if (x(i) <= lower(i) && gradient(i) > 0) free[i] = false;
    if (x(i) >= upper(i) && gradient(i) < 0) free[i] = false;
  }
  return free;
}

double newtonDecrement(const QuadraticModel& model, const std::vector<bool>& free) {
  Eigen::VectorXd g;
  Eigen::MatrixXd H;
  restrict(model, indices(free), g, H);
  if (g.size() == 0) return 0;
  Eigen::LLT<Eigen::MatrixXd> llt(H);
  if (llt.info() != Eigen::Success) return std::numeric_limits<double>::infinity();
  return 0.5 * g.dot(llt.solve(g));
}

bool dampedStep(const QuadraticModel& model, const std::vector<bool>& free, double lambda,
                Eigen::VectorXd& dx) {
  const std::vector<int> idx = indices(free);
  Eigen::VectorXd g;
  Eigen::MatrixXd H;
  restrict(model, idx, g, H);
  dx = Eigen::VectorXd::Zero(model.gradient.size());
  if (g.size() == 0) return true;
  // the floor keeps a parameter with ~zero curvature from taking an undamped step
  const double floor = kDampingFloor * H.diagonal().cwiseAbs().maxCoeff();
  for (int i = 0; i < H.rows(); ++i) H(i, i) += lambda * std::max(std::abs(H(i, i)), floor);
  Eigen::LLT<Eigen::MatrixXd> llt(H);
  if (llt.info() != Eigen::Success) return false;
  const Eigen::VectorXd step = llt.solve(-g);
  if (!step.allFinite()) return false;
  for (size_t i = 0; i < idx.size(); ++i) dx(idx[i]) = step(i);
  return true;
}

double predictedDecrease(const QuadraticModel& model, const Eigen::VectorXd& dx) {
  return -(model.gradient.dot(dx) + 0.5 * dx.dot(model.hessian * dx));
}

LMResult LevenbergMarquardt::minimise(Problem& problem, const Eigen::VectorXd& x0,
                                      const std::function<void(const LMIteration&)>& log) const {
  const Eigen::VectorXd lower = problem.lowerBounds();
  const Eigen::VectorXd upper = problem.upperBounds();

  LMResult result;
  result.x = clampToBox(x0, lower, upper);
  result.chi2 = problem.chi2(result.x);
  result.chi2Calls = 1;
  result.model = problem.model(result.x);
  double lambda = _options.lambda0;
  double lastDecrease = std::numeric_limits<double>::infinity();  ///< chi2 gained by the last step
  std::vector<bool> free(result.x.size(), true);

  for (int iteration = 1; iteration <= _options.maxIterations; ++iteration) {
    free = freeParameters(result.x, result.model.gradient, lower, upper);
    double decrement = newtonDecrement(result.model, free);
    if ((decrement < _options.tolerance || lastDecrease < _options.tolerance) && problem.modelIsApproximate()) {
      // confirm the stopping test with exact derivatives at the same point (README D12), also
      // when the last step gained less than the tolerance; the fit goes on with them
      problem.requestExactModel();
      result.model = problem.model(result.x);
      result.exactModelAt = iteration;
      free = freeParameters(result.x, result.model.gradient, lower, upper);
      decrement = newtonDecrement(result.model, free);
    }
    if (decrement < _options.tolerance) {
      result.converged = true;
      result.message = "Newton decrement below tolerance";
      break;
    }

    bool accepted = false;
    int trials = 0;
    for (int trial = 0; trial < _options.maxTrials && !accepted; ++trial) {
      Eigen::VectorXd dx;
      if (!dampedStep(result.model, free, lambda, dx)) {
        lambda *= 10;
        continue;
      }
      const Eigen::VectorXd xTrial = clampToBox(result.x + dx, lower, upper);
      dx = xTrial - result.x;  // the step actually taken, Eq. (12)
      const double predicted = predictedDecrease(result.model, dx);
      if (!(predicted > 0)) {
        lambda *= 4;
        continue;
      }
      const double chi2Trial = problem.chi2(xTrial);
      ++result.chi2Calls;
      ++trials;
      const double rho = (result.chi2 - chi2Trial) / predicted;
      if (rho > kAcceptRatio) {
        accepted = true;
        lastDecrease = result.chi2 - chi2Trial;
        result.x = xTrial;
        result.chi2 = chi2Trial;
        result.model = problem.model(result.x);  // right after chi2(xTrial), see Problem.h
        result.iterations = iteration;
        if (log) log({iteration, chi2Trial, decrement, predicted, rho, lambda, trials, result.x});
        if (rho > 0.75) lambda = std::max(lambda / 3, kMinLambda);
        if (rho < 0.25) lambda *= 2;
      } else {
        lambda *= 4;
      }
    }
    if (!accepted) {
      result.message = "no acceptable step after " + std::to_string(_options.maxTrials) + " trials";
      break;
    }
    if (iteration == _options.maxIterations) result.message = "maximum number of iterations reached";
  }

  result.atBound.assign(result.x.size(), false);
  for (int i = 0; i < result.x.size(); ++i)
    result.atBound[i] = result.x(i) <= lower(i) || result.x(i) >= upper(i);
  return result;
}

}
}
