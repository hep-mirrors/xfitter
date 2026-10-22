#pragma once

// Damped Gauss-Newton (Levenberg-Marquardt) loop with projected box bounds,
// README.md Eqs. (11)-(12). Pure numerics: no xFitter headers in core/.

#include "core/Problem.h"

#include <functional>
#include <string>
#include <vector>

namespace xfitter {
namespace gn {

struct LMOptions {
  int maxIterations = 50;
  double tolerance = 1e-4;  ///< stop when the Newton decrement, Eq. (11), is below this (chi2 units)
  double lambda0 = 1e-3;
  int maxTrials = 20;       ///< damped steps tried per iteration
};

/// One line of the iteration log
struct LMIteration {
  int iteration;
  double chi2;       ///< after the iteration
  double decrement;  ///< Newton decrement before the step
  double predicted;  ///< predicted decrease of the accepted step
  double rho;        ///< actual / predicted decrease of the accepted step
  double lambda;     ///< damping used for the accepted step
  int trials;        ///< chi2 evaluations in this iteration
  Eigen::VectorXd x; ///< parameters after the iteration
};

struct LMResult {
  Eigen::VectorXd x;
  double chi2 = 0;
  QuadraticModel model;       ///< at x
  int exactModelAt = 0;       ///< iteration at which exact derivatives were requested (0: never)
  std::vector<bool> atBound;  ///< parameter held at a bound by the active set
  int iterations = 0;
  int chi2Calls = 0;
  bool converged = false;
  std::string message;
};

class LevenbergMarquardt {
public:
  explicit LevenbergMarquardt(const LMOptions& options) : _options(options) {}

  LMResult minimise(Problem& problem, const Eigen::VectorXd& x0,
                    const std::function<void(const LMIteration&)>& log = nullptr) const;

private:
  LMOptions _options;
};

/// Parameters not blocked by a bound: a parameter at a bound whose gradient points out of
/// the box is held fixed (active set, Eq. 12).
std::vector<bool> freeParameters(const Eigen::VectorXd& x, const Eigen::VectorXd& gradient,
                                 const Eigen::VectorXd& lower, const Eigen::VectorXd& upper);

/// Newton decrement 1/2 g.H^-1.g over the free parameters, Eq. (11);
/// +infinity if H is not positive definite there.
double newtonDecrement(const QuadraticModel& model, const std::vector<bool>& free);

/// Damped step (H + lambda diag|H|) dx = -g over the free parameters, Eq. (11), with
/// diag|H| floored at 1e-8 max|H_ii|; returns false if the damped matrix is not positive
/// definite or the step is not finite.
bool dampedStep(const QuadraticModel& model, const std::vector<bool>& free, double lambda,
                Eigen::VectorXd& dx);

/// Predicted decrease -(g.dx + 1/2 dx.H.dx), Eq. (11)
double predictedDecrease(const QuadraticModel& model, const Eigen::VectorXd& dx);

}
}
