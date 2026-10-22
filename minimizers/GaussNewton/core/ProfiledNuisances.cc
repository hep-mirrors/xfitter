#include "core/ProfiledNuisances.h"

namespace xfitter {
namespace gn {

ProfiledModel profileNuisances(const PointModel& pm, const Eigen::MatrixXd& jacobianT,
                               const Eigen::VectorXd& priorScale) {
  const Eigen::MatrixXd& B = pm.dResidualdB;                        // m x n, Eq. (7)
  const Eigen::MatrixXd Jp = pm.dResidualdT.asDiagonal() * jacobianT;  // n x P, Eq. (7)

  // Eq. (8); the prior residuals are sqrt(pi_k) b_k
  const Eigen::VectorXd gp = 2 * Jp.transpose() * pm.residual + jacobianT.transpose() * pm.dLogTermdT;
  const Eigen::VectorXd gb = 2 * B * pm.residual + 2 * priorScale.cwiseSqrt().cwiseProduct(pm.prior);
  Eigen::MatrixXd Hpp = 2 * Jp.transpose() * Jp;
  Hpp.noalias() += jacobianT.transpose() * pm.d2LogTermdT2.asDiagonal() * jacobianT;
  const Eigen::MatrixXd Hpb = 2 * Jp.transpose() * B.transpose();
  Eigen::MatrixXd Hbb = 2 * B * B.transpose();
  Hbb.diagonal() += 2 * priorScale;

  // Eq. (9)
  ProfiledModel result;
  result.shiftGradient = gb;
  if (Hbb.rows() == 0) {
    result.model = {gp, Hpp};
    result.shiftCorrection = gb;
    return result;
  }
  const Eigen::LDLT<Eigen::MatrixXd> ldlt(Hbb);
  const Eigen::MatrixXd X = ldlt.solve(Hpb.transpose());  // H_bb^-1 H_bp
  const Eigen::VectorXd y = ldlt.solve(gb);              // H_bb^-1 g_b
  Eigen::MatrixXd H = Hpp - Hpb * X;
  result.model.hessian = 0.5 * (H + H.transpose());
  result.model.gradient = gp - Hpb * y;
  result.shiftCorrection = -y;
  return result;
}

bool covarianceFromHessian(const Eigen::MatrixXd& hessian, Eigen::MatrixXd& covariance) {
  const Eigen::LLT<Eigen::MatrixXd> llt(hessian);
  if (llt.info() == Eigen::Success) {
    covariance = 2 * llt.solve(Eigen::MatrixXd::Identity(hessian.rows(), hessian.cols()));
    return true;
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(hessian);
  const Eigen::VectorXd& values = eigen.eigenvalues();
  Eigen::VectorXd inverse = Eigen::VectorXd::Zero(values.size());
  const double cut = 1e-12 * values.cwiseAbs().maxCoeff();
  for (int i = 0; i < values.size(); ++i)
    if (values(i) > cut) inverse(i) = 1 / values(i);
  covariance = 2 * eigen.eigenvectors() * inverse.asDiagonal() * eigen.eigenvectors().transpose();
  return false;
}

}
}
