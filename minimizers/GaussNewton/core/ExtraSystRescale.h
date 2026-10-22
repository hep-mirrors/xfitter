#pragma once

// ExtraSystRescale, README.md Eq. (13b). xFitter profiles the shifts twice: b1 with the
// nominal errors, then b with the Poisson stat errors rescaled by E_i = exp(-sum_k beta_ki b1_k)
// (linear nuisance sources). The chi2 depends on t also through b1(t). This adds that term to
// the theory gradient and the rescaling to d sigma^2/dt. Pure numerics: no xFitter headers.

#include "core/DataSnapshot.h"

namespace xfitter {
namespace gn {

struct RescaleModel {
  DataSnapshot corrected;     ///< sigma2DT with the rescaled stat part, for computePointModel
  Eigen::VectorXd gradient;   ///< d chi2/dt through the first-pass shifts, added to J_t^T gamma
  Eigen::VectorXd firstShifts;  ///< b1
  double sigmaMismatch = 0;   ///< max |sigma rebuilt from b1 - sigma of xFitter| / sigma
};

/// s: snapshot with extraSystRescale filled
RescaleModel extraSystRescale(const DataSnapshot& s);

}
}
