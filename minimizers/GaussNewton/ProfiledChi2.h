#pragma once

// The xFitter fit as a gn::Problem (README.md, Sec. 2.3):
//   chi2(p)  = official chi2 (fcn + parameter priors), then snapshot of the chi2 state;
//   model(p) = theory Jacobian -> point model -> Schur complement over the shifts, Eqs. (7)-(9).

#include "core/DataSnapshot.h"
#include "core/ProfiledNuisances.h"

namespace xfitter {
namespace gn {
class FortranBridge;
class TheoryEvaluator;
class TheoryJacobian;

/// Consistency checks and timings, accumulated over the fit
struct ProfiledChi2Diagnostics {
  double maxChi2Mismatch = 0;  ///< |module chi2 - official chi2|, Eq. (1)
  double maxShiftChange = 0;   ///< |H_bb^-1 g_b|, ~0 when b is the exact profile, Eq. (9)
  double fcnTime = 0, jacobianTime = 0, algebraTime = 0;
};

class ProfiledChi2 : public Problem {
public:
  ProfiledChi2(TheoryEvaluator& evaluator, const FortranBridge& bridge, TheoryJacobian& jacobian,
               const Eigen::VectorXd& lower, const Eigen::VectorXd& upper);

  int size() const override { return static_cast<int>(_lower.size()); }
  double chi2(const Eigen::VectorXd& p) override;
  QuadraticModel model(const Eigen::VectorXd& p) override;
  Eigen::VectorXd lowerBounds() const override { return _lower; }
  Eigen::VectorXd upperBounds() const override { return _upper; }

  /// fcn flag of the next chi2 call: 1 for the very first call, then 2
  void setFcnFlag(int iflag) { _iflag = iflag; }

  const ProfiledModel& lastProfile() const { return _profile; }
  const ProfiledChi2Diagnostics& diagnostics() const { return _diagnostics; }

private:
  TheoryEvaluator& _evaluator;
  const FortranBridge& _bridge;
  TheoryJacobian& _jacobian;
  Eigen::VectorXd _lower, _upper;
  int _iflag = 2;

  // state of the last chi2 call, see Problem::model
  Eigen::VectorXd _lastP;
  double _lastChi2 = 0;
  DataSnapshot _snapshot;

  ProfiledModel _profile;
  ProfiledChi2Diagnostics _diagnostics;
};

}
}
