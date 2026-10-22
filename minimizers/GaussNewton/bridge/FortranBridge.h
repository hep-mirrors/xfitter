#pragma once

// C++ side of bridge/xfitter_gn_bridge.f: the only place where the module reads the
// xFitter chi2 state. Read-only. See README.md, Sec. 2.

#include "core/DataSnapshot.h"

#include <vector>

namespace xfitter {
namespace gn {

struct ChiSquareConfig {
  int fitPoints = 0;          ///< fitted points with diagonal errors
  int nuisances = 0;          ///< isNuisance sources
  int covariancePoints = 0;   ///< points treated with a covariance matrix
  int matrixSources = 0;
  int offsetSources = 0;
  int externalSources = 0;
  bool poissonCorr = false;
  bool extraSystRescale = false;
  bool firstIterationRescale = false;
  bool asymmetricErrors = false;
};

class FortranBridge {
public:
  /// Reads the configuration; call after the data are read (atStart or later)
  FortranBridge();

  const ChiSquareConfig& config() const { return _config; }

  /// Stops on configurations the module cannot handle and warns on those where the
  /// derivative model is approximate (README.md, Sec. 3, Eq. 13)
  void checkSupported() const;

  /// Chi2 ingredients of the last fcn call. Take it right after fcn, before any
  /// theory-only evaluation changes THEO.
  DataSnapshot snapshot() const;

  /// THEO of the fitted points, e.g. after update_theory_iteration
  Eigen::VectorXd theory() const;

  /// Polarisation-shift terms that fcn adds to the chi2
  double polarisationChi2() const;

private:
  ChiSquareConfig _config;
  std::vector<int> _points;   ///< 1-based xFitter indices of the fitted points
  std::vector<int> _sources;  ///< 1-based xFitter indices of the nuisance sources
};

}
}
