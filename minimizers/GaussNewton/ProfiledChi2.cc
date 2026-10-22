#include "ProfiledChi2.h"
#include "bridge/FortranBridge.h"
#include "bridge/TheoryEvaluator.h"
#include "core/PointModel.h"
#include "jacobian/TheoryJacobian.h"
#include "util/Timer.h"
#include "xfitter_cpp_base.h"

#include <algorithm>
#include <cmath>

namespace xfitter {
namespace gn {

ProfiledChi2::ProfiledChi2(TheoryEvaluator& evaluator, const FortranBridge& bridge,
                           TheoryJacobian& jacobian, const Eigen::VectorXd& lower,
                           const Eigen::VectorXd& upper)
    : _evaluator(evaluator), _bridge(bridge), _jacobian(jacobian), _lower(lower), _upper(upper) {}

double ProfiledChi2::chi2(const Eigen::VectorXd& p) {
  Timer::Scope timer(_diagnostics.fcnTime);
  _lastChi2 = _evaluator.chi2(p, _iflag);
  _iflag = 2;
  // the snapshot must be taken now: the Jacobian changes THEO
  _snapshot = _bridge.snapshot();
  _lastP = p;
  return _lastChi2;
}

QuadraticModel ProfiledChi2::model(const Eigen::VectorXd& p) {
  if (_lastP.size() != p.size() || _lastP != p) chi2(p);  // see Problem::model

  Eigen::MatrixXd J;
  {
    Timer::Scope timer(_diagnostics.jacobianTime);
    J = _jacobian.compute(p, _snapshot.theory);
  }
  if (!J.allFinite()) {
    // e.g. a perturbed parameter gives an invalid PDF: these derivatives are dropped
    hf_errlog(2026102231, "W: GaussNewton: non-finite theory derivatives set to zero");
    J = J.unaryExpr([](double v) { return std::isfinite(v) ? v : 0.; });
  }

  Timer::Scope timer(_diagnostics.algebraTime);
  const PointModel pm = computePointModel(_snapshot);
  const double mismatch = std::abs(pm.chi2() + _evaluator.priorChi2(p) - _lastChi2);
  _diagnostics.maxChi2Mismatch = std::max(_diagnostics.maxChi2Mismatch, mismatch);
  if (mismatch > 1e-6 * std::max(1., _lastChi2))
    hf_errlog(2026102230, "W: GaussNewton: module chi2 differs from fcn, the derivative model may be wrong");

  _profile = profileNuisances(pm, J, _snapshot.priorScale);
  const double shiftChange = _profile.shiftCorrection.size() ? _profile.shiftCorrection.cwiseAbs().maxCoeff() : 0;
  _diagnostics.maxShiftChange = std::max(_diagnostics.maxShiftChange, shiftChange);

  QuadraticModel model = _profile.model;
  _evaluator.addPriors(p, model);
  return model;
}

}
}
