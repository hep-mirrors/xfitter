#pragma once

// Gauss-Newton model of the chi2 on the joint vector (p, b) and its profile over the
// nuisance shifts b by a Schur complement, README.md Eqs. (7)-(10).
// Pure numerics: no xFitter headers in core/.

#include "core/PointModel.h"
#include "core/Problem.h"

namespace xfitter {
namespace gn {

struct ProfiledModel {
  QuadraticModel model;           ///< g_eff and H_eff over p, Eq. (9)
  Eigen::VectorXd shiftGradient;  ///< g_b: ~0 when b is xFitter's exact profile, Eq. (9)
  Eigen::VectorXd shiftCorrection;  ///< -H_bb^-1 g_b: change of b predicted by the model
};

/// jacobianT: dt_i/dp_j, n points x n parameters
ProfiledModel profileNuisances(const PointModel& pointModel, const Eigen::MatrixXd& jacobianT,
                               const Eigen::VectorXd& priorScale);

/// C = 2 H^-1, Eq. (10). If H is not positive definite, C is the pseudo-inverse over the
/// positive eigenvalues and the function returns false.
bool covarianceFromHessian(const Eigen::MatrixXd& hessian, Eigen::MatrixXd& covariance);

}
}
