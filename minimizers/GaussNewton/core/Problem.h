#pragma once

// Interface between the Levenberg-Marquardt loop and the fit (README.md, Sec. 3).
// Pure numerics: no xFitter headers in core/.

#include <Eigen/Dense>

#include <limits>

namespace xfitter {
namespace gn {

/// Local quadratic model: chi2(x + dx) ~ chi2(x) + g.dx + 1/2 dx.H.dx (chi2 units)
struct QuadraticModel {
  Eigen::VectorXd gradient;
  Eigen::MatrixXd hessian;
};

class Problem {
public:
  virtual ~Problem() = default;

  /// Number of parameters
  virtual int size() const = 0;

  /// The objective at x. For xFitter this is always the official chi2 from fcn.
  virtual double chi2(const Eigen::VectorXd& x) = 0;

  /// Quadratic model at x. The loop calls it only right after chi2(x) at the same x,
  /// so an implementation may reuse the state of that evaluation.
  virtual QuadraticModel model(const Eigen::VectorXd& x) = 0;

  /// Box bounds, +-infinity where there is none
  virtual Eigen::VectorXd lowerBounds() const {
    return Eigen::VectorXd::Constant(size(), -std::numeric_limits<double>::infinity());
  }
  virtual Eigen::VectorXd upperBounds() const {
    return Eigen::VectorXd::Constant(size(), std::numeric_limits<double>::infinity());
  }
};

}
}
