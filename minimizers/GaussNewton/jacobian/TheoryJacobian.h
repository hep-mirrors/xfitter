#pragma once

// Interface for J_t = dt/dp, the Jacobian of the theory of the fitted points
// (README.md, Eq. 7). Implementations: FiniteDifferenceJacobian, FactorisedJacobian.

#include <Eigen/Dense>

namespace xfitter {
namespace gn {

class TheoryJacobian {
public:
  virtual ~TheoryJacobian() = default;

  /// J_t at p (n points x n parameters), given t0 = t(p). May change the global xFitter
  /// state: call it after the snapshot of the fcn call at p has been taken.
  virtual Eigen::MatrixXd compute(const Eigen::VectorXd& p, const Eigen::VectorXd& t0) = 0;

  /// Theory evaluations used so far
  virtual int evaluations() const = 0;

  /// True if the last Jacobian is an approximation
  virtual bool isApproximate() const { return false; }

  /// From now on compute() gives exact (finite-difference) Jacobians
  virtual void requestExact() {}

  /// Called after every chi2 evaluation with its point and theory t(p)
  virtual void observe(const Eigen::VectorXd& p, const Eigen::VectorXd& t) {}
};

}
}
