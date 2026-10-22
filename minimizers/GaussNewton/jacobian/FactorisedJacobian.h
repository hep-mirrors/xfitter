#pragma once

// J_t = K F, README.md Eqs. (14)-(17): F = df0/dp of the input PDFs (cheap, no
// evolution), K = dt/df0 known on the input-PDF directions explored so far. Columns whose
// F_j is not in the explored span, and the columns of parameters that do not change the
// input PDFs, come from finite differences; every theory evaluation adds an exact pair.

#include "core/LinearResponse.h"
#include "jacobian/FiniteDifferenceJacobian.h"

#include <vector>

namespace xfitter {
namespace gn {
class InputPdfs;

struct FactorisedOptions {
  double spanTolerance = 1e-3;    ///< column j from K F_j if |F_j - Q Q^T F_j| <= this |F_j|, Eq. (17)
  double addTolerance = 1e-4;     ///< a pair adds a direction if |u - Q Q^T u| > this |u|, Eq. (15)
  double verifyTolerance = 2e-2;  ///< linearity check: |J_KF - J_FD| / |J_FD| (sigma units)
};

struct FactorisedStatistics {
  int fdColumns = 0;            ///< columns from finite differences
  int responseColumns = 0;      ///< columns from K F
  int resets = 0;               ///< K dropped because a non-PDF parameter changed
  double verifyDeviation = -1;  ///< last linearity check, -1: none since the start or the last reset
  bool linear = true;           ///< false after a failed linearity check
};

class FactorisedJacobian : public TheoryJacobian {
public:
  FactorisedJacobian(FiniteDifferenceJacobian& fd, InputPdfs& input, const FactorisedOptions& options);

  /// Point uncertainties, the metric of the linearity check (default: 1)
  void setPointScale(const Eigen::VectorXd& sigma) { _sigma = sigma; }

  Eigen::MatrixXd compute(const Eigen::VectorXd& p, const Eigen::VectorXd& t0) override;
  int evaluations() const override { return _fd.evaluations(); }
  /// The pair (f0(p) - f0(previous), t - t(previous)), Eq. (15)
  void observe(const Eigen::VectorXd& p, const Eigen::VectorXd& t) override;

  const FactorisedStatistics& statistics() const { return _stats; }
  int dimension() const { return _response.dimension(); }

private:
  /// F = df0/dp by central differences (one-sided at a bound), Eq. (14)
  Eigen::MatrixXd inputJacobian(const Eigen::VectorXd& p, const Eigen::VectorXd& f);
  /// Drops K if a parameter outside the input PDFs moved since the last call
  void checkNonPdfParameters(const Eigen::VectorXd& p);
  /// Linearity check, design decision D12: column j by finite differences, compared with K F_j.
  /// Returns the finite-difference column.
  Eigen::VectorXd verify(const Eigen::VectorXd& p, const Eigen::VectorXd& t0, int j, const Eigen::VectorXd& column);
  /// The exact pair of the finite-difference column j at p, Eq. (15)
  void addColumnPair(const Eigen::VectorXd& p, const Eigen::VectorXd& f, int j, const Eigen::VectorXd& column);

  FiniteDifferenceJacobian& _fd;
  InputPdfs& _input;
  Eigen::VectorXd _sigma;
  FactorisedOptions _options;
  LinearResponse _response;
  std::vector<bool> _nonPdf;      ///< parameter does not change f0; set at the first compute
  Eigen::VectorXd _pLast, _fLast, _tLast;  ///< last observed point
  Eigen::VectorXd _pNonPdf;       ///< parameters at which K is valid
  FactorisedStatistics _stats;
};

}
}
