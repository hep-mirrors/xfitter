#include "core/PointModel.h"

#include <algorithm>
#include <cmath>

namespace xfitter {
namespace gn {

double dSigma2dT(double d, double t, double sigma2DT, double sigma2TT) {
  if (t <= 0) return 0;
  return sigma2DT * std::abs(d) + 2 * sigma2TT * t;
}

PointModel computePointModel(const DataSnapshot& s) {
  const int n = s.points();
  const int m = s.sources();
  PointModel pm;
  pm.residual.resize(n);
  pm.dResidualdT.resize(n);
  pm.dResidualdB.resize(m, n);
  pm.logTerm = Eigen::VectorXd::Zero(n);
  pm.dLogTermdT = Eigen::VectorXd::Zero(n);
  pm.d2LogTermdT2 = Eigen::VectorXd::Zero(n);
  pm.prior = s.priorScale.cwiseSqrt().cwiseProduct(s.shifts);

  for (int i = 0; i < n; ++i) {
    const double d = s.data(i), t = s.theory(i), sigma = s.sigma(i);
    const double sigma2 = sigma * sigma;

    // Eq. (4)
    const Eigen::VectorXd gamma = s.gamma.col(i);
    pm.residual(i) = (d - t + gamma.dot(s.shifts)) / sigma;
    pm.dResidualdB.col(i) = gamma / sigma;

    // Eq. (5): dGamma_ki/dt_i = power_k Gamma_ki / t_i
    double dShiftdT = 0;
    if (!s.additive[i] && t != 0)
      dShiftdT = s.gammaPower.cwiseProduct(gamma).dot(s.shifts) / t;
    const double ds2 = dSigma2dT(d, t, s.sigma2DT(i), s.sigma2TT(i));
    pm.dResidualdT(i) = (-1 + dShiftdT) / sigma - pm.residual(i) * ds2 / (2 * sigma2);

    // Eq. (6): L = log(sigma^2/alpha^2), keep only the non-negative curvature
    if (s.logTerm[i]) {
      const double d2s2 = t > 0 ? 2 * s.sigma2TT(i) : 0;
      pm.logTerm(i) = std::log(sigma2 / (s.alpha(i) * s.alpha(i)));
      pm.dLogTermdT(i) = ds2 / sigma2;
      pm.d2LogTermdT2(i) = std::max(0., d2s2 / sigma2 - std::pow(ds2 / sigma2, 2));
    }
  }
  return pm;
}

}
}
