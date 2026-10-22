#include "core/LinearResponse.h"

namespace xfitter {
namespace gn {

void LinearResponse::project(const Eigen::VectorXd& u, Eigen::VectorXd& a, Eigen::VectorXd& r) const {
  r = u;
  a = Eigen::VectorXd::Zero(_Q.cols());
  if (_Q.cols() == 0) return;
  for (int pass = 0; pass < 2; ++pass) {  // twice is enough for orthogonality to round-off
    const Eigen::VectorXd c = _Q.transpose() * r;
    r -= _Q * c;
    a += c;
  }
}

bool LinearResponse::add(const Eigen::VectorXd& u, const Eigen::VectorXd& w) {
  if (!u.allFinite() || !w.allFinite()) return false;
  if (_Q.cols() > 0 && (u.size() != _Q.rows() || w.size() != _Z.rows())) return false;
  const double norm = u.norm();
  if (!(norm > 0)) return false;
  Eigen::VectorXd a, r;
  project(u, a, r);
  const double rn = r.norm();
  if (!(rn > _addTolerance * norm)) return false;
  // K q = (K u - K Q a) / |r| = (w - Z a) / |r|, Eq. (15)
  const Eigen::VectorXd z = _Q.cols() > 0 ? Eigen::VectorXd((w - _Z * a) / rn) : Eigen::VectorXd(w / rn);
  const int m = static_cast<int>(_Q.cols());
  _Q.conservativeResize(u.size(), m + 1);
  _Z.conservativeResize(w.size(), m + 1);
  _Q.col(m) = r / rn;
  _Z.col(m) = z;
  return true;
}

LinearResponse::Projection LinearResponse::apply(const Eigen::VectorXd& u) const {
  Projection result;
  Eigen::VectorXd a, r;
  project(u, a, r);
  result.response = _Q.cols() > 0 ? Eigen::VectorXd(_Z * a) : Eigen::VectorXd();
  const double norm = u.norm();
  result.residual = norm > 0 ? r.norm() / norm : 0.;
  return result;
}

void LinearResponse::clear() {
  _Q.resize(0, 0);
  _Z.resize(0, 0);
}

}
}
