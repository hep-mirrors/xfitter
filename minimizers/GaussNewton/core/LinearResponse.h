#pragma once

// The linear response K = dt/df0 of the theory to the input PDFs, known on a subspace
// (README.md, Eqs. 15-16). Pure numerics: Eigen + std only.
//
// The response is stored as an orthonormal basis Q of the explored input-PDF directions
// and Z = K Q. Each exact pair (u, w = K u), e.g. the input PDFs and theories of two
// evaluated points, adds one direction unless u is already in the span.

#include <Eigen/Dense>

namespace xfitter {
namespace gn {

class LinearResponse {
public:
  /// addTolerance: a pair is added only if |u - Q Q^T u| > addTolerance |u|, Eq. (15)
  explicit LinearResponse(double addTolerance = 1e-6) : _addTolerance(addTolerance) {}

  /// Adds the exact pair (u, w = K u). Returns false if u is in the span or not usable.
  bool add(const Eigen::VectorXd& u, const Eigen::VectorXd& w);

  struct Projection {
    Eigen::VectorXd response;  ///< K Q Q^T u, Eq. (16)
    double residual = 0;       ///< |u - Q Q^T u| / |u| (0 for u = 0)
  };
  /// K applied to the part of u inside the span
  Projection apply(const Eigen::VectorXd& u) const;

  int dimension() const { return static_cast<int>(_Q.cols()); }
  void clear();

private:
  /// coefficients a = Q^T u and residual r = u - Q a, Gram-Schmidt applied twice
  void project(const Eigen::VectorXd& u, Eigen::VectorXd& a, Eigen::VectorXd& r) const;

  double _addTolerance;
  Eigen::MatrixXd _Q;  ///< n_input x m, orthonormal
  Eigen::MatrixXd _Z;  ///< n_points x m, Z = K Q
};

}
}
