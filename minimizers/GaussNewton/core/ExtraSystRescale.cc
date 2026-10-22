#include "core/ExtraSystRescale.h"
#include "core/PointModel.h"

#include <algorithm>
#include <cmath>

namespace xfitter {
namespace gn {

RescaleModel extraSystRescale(const DataSnapshot& s) {
  const int n = s.points(), m = s.sources();
  RescaleModel result;

  // first pass: nominal errors, exact profile b1 (r is linear in b)
  DataSnapshot first = s;
  first.sigma = s.sigmaFirstPass;
  first.shifts = Eigen::VectorXd::Zero(m);
  const PointModel at0 = computePointModel(first);
  const Eigen::MatrixXd& B1 = at0.dResidualdB;  // m x n
  Eigen::MatrixXd A = B1 * B1.transpose();
  A.diagonal() += s.priorScale;
  const Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
  result.firstShifts = m > 0 ? Eigen::VectorXd(-ldlt.solve(B1 * at0.residual)) : Eigen::VectorXd();
  first.shifts = result.firstShifts;
  const PointModel pm1 = computePointModel(first);

  // rescaling of the stat errors and its check against xFitter's sigma
  const Eigen::VectorXd E = m > 0 ? Eigen::VectorXd((-(s.betaLinear.transpose() * result.firstShifts)).array().exp())
                                  : Eigen::VectorXd::Ones(n);
  Eigen::VectorXd stat(n);  // Poisson stat part of sigma^2 before the rescaling
  for (int i = 0; i < n; ++i) {
    const double t = s.theory(i) > 0 ? s.theory(i) : s.data(i);  // as GetPointErrors
    stat(i) = s.sigma2DTStat(i) * std::abs(s.data(i) * t);
    const double rebuilt = std::sqrt(s.sigmaFirstPass(i) * s.sigmaFirstPass(i) + stat(i) * (E(i) - 1));
    result.sigmaMismatch = std::max(result.sigmaMismatch, std::abs(rebuilt - s.sigma(i)) / s.sigma(i));
  }
  result.corrected = s;
  result.corrected.sigma2DT = s.sigma2DT + s.sigma2DTStat.cwiseProduct(E - Eigen::VectorXd::Ones(n));

  // d chi2/d sigma_j^2 at the final shifts, and d sigma_j^2 / d b1_k = -beta_kj stat_j E_j
  const PointModel pm2 = computePointModel(result.corrected);
  Eigen::VectorXd c(n);
  for (int i = 0; i < n; ++i) {
    const double s2 = s.sigma(i) * s.sigma(i);
    c(i) = (-pm2.residual(i) * pm2.residual(i) + (s.logTerm[i] ? 1. : 0.)) / s2;
  }
  const Eigen::VectorXd w = -s.betaLinear * c.cwiseProduct(stat).cwiseProduct(E);  // m

  // d b1/dt = -Phi1_bb^-1 Phi1_bt with Phi1_bb = 2 A and, per point,
  // Phi1_bt(k, i) = 2 B1_ki a1_i + 2 r1_i d2r1_i/dt_i db_k
  result.gradient = Eigen::VectorXd::Zero(n);
  if (m == 0) return result;
  const Eigen::VectorXd y = ldlt.solve(w) / 2;
  for (int i = 0; i < n; ++i) {
    const double t = s.theory(i), s1 = s.sigmaFirstPass(i) * s.sigmaFirstPass(i);
    const double ds1 = dSigma2dT(s.data(i), t, s.sigma2DT(i), s.sigma2TT(i));
    double sum = 0;
    for (int k = 0; k < m; ++k) {
      const double dGammadT = (!s.additive[i] && t != 0) ? s.gammaPower(k) * s.gamma(k, i) / t : 0.;
      const double d2r = dGammadT / s.sigmaFirstPass(i) - B1(k, i) * ds1 / (2 * s1);
      sum += y(k) * (2 * B1(k, i) * pm1.dResidualdT(i) + 2 * pm1.residual(i) * d2r);
    }
    result.gradient(i) = -sum;
  }
  return result;
}

}
}
