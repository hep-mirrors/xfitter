#pragma once

// The xFitter chi2 ingredients at one parameter point, as plain arrays (README.md, Sec. 3).
// Filled by bridge/FortranBridge right after an fcn call. Only fitted, diagonal points and
// nuisance-parameter sources are included. Pure data: no xFitter headers in core/.

#include <Eigen/Dense>

#include <vector>

namespace xfitter {
namespace gn {

struct DataSnapshot {
  // per point i (n points)
  Eigen::VectorXd data;          ///< d_i
  Eigen::VectorXd theory;        ///< t_i
  Eigen::VectorXd sigma;         ///< sigma_i(t) as used by xFitter, Eq. (3)
  Eigen::VectorXd sigma2DT;      ///< coefficient of |d_i t_i| in sigma_i^2, Eq. (3)
  Eigen::VectorXd sigma2TT;      ///< coefficient of t_i^2 in sigma_i^2, Eq. (3)
  std::vector<bool> additive;    ///< all sources additive for this point (LForceAdditiveData)
  std::vector<bool> logTerm;     ///< point enters the log term, Eq. (6)
  Eigen::VectorXd alpha;         ///< alpha_i of the log term, Eq. (6)

  // per source k (m sources) and point
  Eigen::MatrixXd gamma;         ///< Gamma_ki, m x n, Eq. (2)
  Eigen::VectorXd gammaPower;    ///< exponent of t in Gamma_ki: 0, 1/2 or 1, Eq. (2)
  Eigen::VectorXd shifts;        ///< b_k, profiled by xFitter
  Eigen::VectorXd priorScale;    ///< pi_k, Eq. (1)

  int points() const { return static_cast<int>(data.size()); }
  int sources() const { return static_cast<int>(shifts.size()); }
};

}
}
