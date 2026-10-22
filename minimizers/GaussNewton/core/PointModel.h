#pragma once

// Per-point residuals and their derivatives with respect to the theory and the
// nuisance shifts, README.md Eqs. (2)-(6). Pure numerics: no xFitter headers in core/.

#include "core/DataSnapshot.h"

namespace xfitter {
namespace gn {

struct PointModel {
  Eigen::VectorXd residual;      ///< r_i, Eq. (4)
  Eigen::VectorXd dResidualdT;   ///< a_i = dr_i/dt_i, Eq. (5)
  Eigen::MatrixXd dResidualdB;   ///< Gamma_ki / sigma_i, m x n, Eq. (5)
  Eigen::VectorXd logTerm;       ///< L_i, zero where the point has no log term, Eq. (6)
  Eigen::VectorXd dLogTermdT;    ///< dL_i/dt_i, Eq. (6)
  Eigen::VectorXd d2LogTermdT2;  ///< max(d2L_i/dt_i^2, 0), Eq. (6)
  Eigen::VectorXd prior;         ///< sqrt(pi_k) b_k, the prior residuals, Eq. (1)

  /// Sum r^2 + sum pi b^2 + sum L, Eq. (1) without parameter priors
  double chi2() const { return residual.squaredNorm() + prior.squaredNorm() + logTerm.sum(); }
};

PointModel computePointModel(const DataSnapshot& snapshot);

/// d sigma_i^2 / d t_i, Eq. (3); zero where xFitter replaces t by d (t <= 0)
double dSigma2dT(double d, double t, double sigma2DT, double sigma2TT);

}
}
