#pragma once

// Observed Hessian of the objective by second differences, as MINUIT HESSE does
// (README.md, Eq. 10 and design decision D10). Pure numerics: no xFitter headers in core/.

#include "core/Problem.h"

namespace xfitter {
namespace gn {

struct ObservedHessianResult {
  Eigen::MatrixXd hessian;
  int chi2Calls = 0;
};

/// Second differences of problem.chi2 around x, where chi2(x) = chi2x.
/// Steps: h_j such that the model curvature predicts a change targetDelta of chi2,
/// refined once from the measured diagonal curvature. Diagonal: central differences
/// (one-sided at a bound). Off-diagonal: forward, f(x+h_j+h_k) with the diagonal points,
/// error O(h), P(P+3)/2 calls in total, as HESSE; or central, also f(x-h_j-h_k), error
/// O(h^2), P(P+1) calls (forward for a pair with a parameter at a bound).
ObservedHessianResult observedHessian(Problem& problem, const Eigen::VectorXd& x, double chi2x,
                                      const Eigen::MatrixXd& modelHessian, double targetDelta = 1e-3,
                                      bool centralOffDiagonal = false);

}
}
