#pragma once

// checkDerivatives: compares the model gradient and Hessian with finite differences of the
// official chi2 at one point (README.md, Sec. 5). This tests Eqs. (2)-(9) against fcn.

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace xfitter {
namespace gn {
class ProfiledChi2;

/// Prints the comparison; returns the largest relative deviation of the gradient
double checkDerivatives(ProfiledChi2& problem, const Eigen::VectorXd& p, const Eigen::VectorXd& steps,
                        const std::vector<std::string>& names);

}
}
