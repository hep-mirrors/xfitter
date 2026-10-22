#include "jacobian/FactorisedJacobian.h"
#include "bridge/InputPdfs.h"
#include "xfitter_cpp_base.h"

#include <cmath>
#include <cstdio>

namespace xfitter {
namespace gn {

FactorisedJacobian::FactorisedJacobian(FiniteDifferenceJacobian& fd, InputPdfs& input,
                                       const FactorisedOptions& options)
    : _fd(fd), _input(input), _options(options) {}

Eigen::MatrixXd FactorisedJacobian::inputJacobian(const Eigen::VectorXd& p, const Eigen::VectorXd& f) {
  Eigen::MatrixXd F(f.size(), p.size());
  for (int j = 0; j < p.size(); ++j) {
    const FiniteDifferenceJacobian::Displacement d = _fd.displacement(p, j, true);
    if (d.up == d.down) {
      F.col(j).setZero();
      continue;
    }
    Eigen::VectorXd q = p;
    q(j) = p(j) + d.up;
    const Eigen::VectorXd up = _input.sample(q);
    q(j) = p(j) + d.down;
    F.col(j) = (up - (d.down == 0 ? f : _input.sample(q))) / (d.up - d.down);
  }
  return F;
}

void FactorisedJacobian::checkNonPdfParameters(const Eigen::VectorXd& p) {
  bool moved = false;
  if (_pNonPdf.size() == p.size())
    for (int j = 0; j < p.size(); ++j) moved = moved || (_nonPdf[j] && p(j) != _pNonPdf(j));
  if (moved) {
    // K depends on alpha_s, masses, normalisations, ...: the stored response is stale
    _response.clear();
    _pLast.resize(0);
    ++_stats.resets;
  }
  _pNonPdf = p;
}

void FactorisedJacobian::observe(const Eigen::VectorXd& p, const Eigen::VectorXd& t) {
  if (!_stats.linear || !t.allFinite()) return;
  const Eigen::VectorXd f = _input.sample(p);
  if (!f.allFinite()) return;
  if (_nonPdf.empty()) {
    const Eigen::MatrixXd F = inputJacobian(p, f);
    for (int j = 0; j < p.size(); ++j) _nonPdf.push_back(F.col(j).cwiseAbs().maxCoeff() == 0);
  }
  checkNonPdfParameters(p);
  if (_pLast.size() == p.size()) _response.add(f - _fLast, t - _tLast);
  _pLast = p;
  _fLast = f;
  _tLast = t;
}

void FactorisedJacobian::verify(const Eigen::VectorXd& p, const Eigen::VectorXd& t0, int j,
                                const Eigen::VectorXd& column) {
  const Eigen::VectorXd exact = _fd.columns(p, t0, {j}).col(0);
  const Eigen::VectorXd scale = _sigma.size() == t0.size() ? _sigma : Eigen::VectorXd::Ones(t0.size());
  const double norm = exact.cwiseQuotient(scale).norm();
  _stats.verifyDeviation = norm > 0 ? (column - exact).cwiseQuotient(scale).norm() / norm : 0.;
  std::printf("Factorised Jacobian: linearity check on parameter %d: |J_KF - J_FD| / |J_FD| = %.2e\n", j,
              _stats.verifyDeviation);
  if (!(_stats.verifyDeviation <= _options.verifyTolerance)) {
    hf_errlog(2026102270, "W: GaussNewton: the theory is not linear in the input PDFs "
                          "(e.g. hadronic data), factorised Jacobian replaced by finite differences");
    _stats.linear = false;
    _response.clear();
  }
}

Eigen::MatrixXd FactorisedJacobian::compute(const Eigen::VectorXd& p, const Eigen::VectorXd& t0) {
  const int P = static_cast<int>(p.size());
  if (_pLast.size() != P || _pLast != p) observe(p, t0);
  if (!_stats.linear || _pLast.size() != P || _pLast != p) {  // no usable f0(p)
    _stats.fdColumns += P;
    return _fd.compute(p, t0);
  }
  const Eigen::VectorXd f = _fLast;
  const Eigen::MatrixXd F = inputJacobian(p, f);

  // Eqs. (16)-(17): K F_j where F_j is in the explored span
  Eigen::MatrixXd J(t0.size(), P);
  std::vector<int> fdColumns;
  int check = -1;
  double checkNorm = -1;
  for (int j = 0; j < P; ++j) {
    LinearResponse::Projection projection;
    if (!_nonPdf[j] && _response.dimension() > 0) projection = _response.apply(F.col(j));
    // also a non-finite F_j (invalid PDF at a displaced point) gives a finite-difference column
    if (_nonPdf[j] || _response.dimension() == 0 || !(projection.residual <= _options.spanTolerance)) {
      fdColumns.push_back(j);
      continue;
    }
    J.col(j) = projection.response;
    const double norm = _sigma.size() == t0.size() ? projection.response.cwiseQuotient(_sigma).norm()
                                                   : projection.response.norm();
    if (norm > checkNorm) {
      check = j;
      checkNorm = norm;
    }
  }
  if (_stats.verifyDeviation < 0 && check >= 0) {
    verify(p, t0, check, J.col(check));
    if (!_stats.linear) {
      _stats.fdColumns += P;
      return _fd.compute(p, t0);
    }
  }
  _stats.responseColumns += P - static_cast<int>(fdColumns.size());
  _stats.fdColumns += static_cast<int>(fdColumns.size());
  if (fdColumns.empty()) return J;

  // the remaining columns by finite differences; each gives an exact pair, Eq. (15)
  const Eigen::MatrixXd columns = _fd.columns(p, t0, fdColumns);
  for (size_t k = 0; k < fdColumns.size(); ++k) {
    const int j = fdColumns[k];
    J.col(j) = columns.col(k);
    const FiniteDifferenceJacobian::Displacement d = _fd.displacement(p, j);
    if (_nonPdf[j] || d.up == d.down) continue;
    Eigen::VectorXd q = p;
    q(j) = p(j) + d.up;
    const Eigen::VectorXd up = _input.sample(q);
    q(j) = p(j) + d.down;
    _response.add(up - (d.down == 0 ? f : _input.sample(q)), columns.col(k) * (d.up - d.down));
  }
  return J;
}

}
}
