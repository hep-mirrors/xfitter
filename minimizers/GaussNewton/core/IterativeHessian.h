#pragma once

// Iterative Hessian of Pumplin et al., hep-ph/0008191, as MINUIT ITERATE (README.md,
// Eq. 18): the Hessian is measured again in its own eigenvector coordinates, scaled so that
// the chi2 rises by delta along each of them, until the measurement is the identity.
// Pure numerics: no xFitter headers in core/.

#include "core/Problem.h"

namespace xfitter {
namespace gn {

struct IterativeHessianOptions {
  double delta = 1;       ///< chi2 rise of the steps along the eigenvectors (MINUIT UP)
  int maxIterations = 10;
  double tolerance = 0.02;  ///< stop when max |H_y / 2 - 1| < tolerance in the eigenvector coordinates
  bool central = true;      ///< central off-diagonal differences
  int maxCycles = 6;        ///< batches to find, per eigenvector, the step with Delta chi2 = delta
};

struct IterativeHessianResult {
  Eigen::MatrixXd hessian;
  /// Eigenvector displacements for Delta chi2 = 1 in the quadratic approximation, one per
  /// column, largest curvature first (MINUIT umat)
  Eigen::MatrixXd displacements;
  int iterations = 0;
  int chi2Calls = 0;
  double deviation = 0;  ///< max |H_y / 2 - 1| of the iteration that gave the result
  int droppedPairs = 0;  ///< off-diagonal pairs with an invalid or far point (correlation not measured)
  bool converged = false;
  double lowestChi2 = 0;  ///< lowest chi2 seen; below the start value: x was not the minimum
};

/// startHessian: e.g. the Gauss-Newton or observed Hessian at x, where chi2(x) = chi2x.
/// The steps ignore the bounds of the problem. Without convergence the result is the
/// Hessian of the most self-consistent iteration.
IterativeHessianResult iterativeHessian(Problem& problem, const Eigen::VectorXd& x, double chi2x,
                                        const Eigen::MatrixXd& startHessian,
                                        const IterativeHessianOptions& options);

/// Displacements x -> x + U_j for Delta chi2 = 1 of a positive definite H (chi2 = 1/2 dx H dx):
/// U_j = v_j sqrt(2 / epsilon_j), largest epsilon first, largest component of v_j positive.
/// Non-positive eigenvalues are replaced by 1e-12 max epsilon. Returns false if there was one.
bool eigenDisplacements(const Eigen::MatrixXd& hessian, Eigen::MatrixXd& displacements);

}
}
