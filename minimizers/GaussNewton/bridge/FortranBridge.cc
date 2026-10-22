#include "bridge/FortranBridge.h"
#include "xfitter_cpp_base.h"

extern "C" {
  void gn_bridge_sizes_(int& nFit, int& nNuis, int& nCov, int& nMatrix, int& nOffset,
                        int& nExtern, int* flags);
  void gn_bridge_indices_(const int& nFit, int* iPoint, const int& nNuis, int* iSource);
  void gn_bridge_theory_(const int& nFit, const int* iPoint, double* t);
  void gn_bridge_snapshot_(const int& nFit, const int* iPoint, const int& nNuis, const int* iSource,
                           double* d, double* t, double* sigma, double* s2dt, double* s2tt,
                           int* additive, int* logTerm, double* alpha, double* gamma,
                           double* power, double* shift, double* prior);
  void gn_bridge_rescale_(const int& nFit, const int* iPoint, const int& nNuis, const int* iSource,
                          double* sigma1, double* s2dtStat, double* betaLinear);
  void gn_bridge_polarisation_chi2_(double& chi2);
}

namespace xfitter {
namespace gn {

FortranBridge::FortranBridge() {
  int flags[4];
  ChiSquareConfig& c = _config;
  gn_bridge_sizes_(c.fitPoints, c.nuisances, c.covariancePoints, c.matrixSources,
                   c.offsetSources, c.externalSources, flags);
  c.poissonCorr = flags[0];
  c.extraSystRescale = flags[1];
  c.firstIterationRescale = flags[2];
  c.asymmetricErrors = flags[3];
  _points.resize(c.fitPoints);
  _sources.resize(c.nuisances);
  gn_bridge_indices_(c.fitPoints, _points.data(), c.nuisances, _sources.data());
}

void FortranBridge::checkSupported() const {
  const ChiSquareConfig& c = _config;
  if (c.covariancePoints > 0)
    hf_errlog(2026102211, "F: GaussNewton does not support covariance-matrix points yet, "
                          "use CovarToNuisance");
  // External sources are fit parameters that enter the chi2 only through the shifts, not the
  // theory: the derivative model has no gradient for them
  if (c.externalSources > 0)
    hf_errlog(2026102214, "F: GaussNewton does not support External systematic sources");
  if (c.offsetSources > 0)
    hf_errlog(2026102212, "W: GaussNewton: Offset sources are not in the derivative model, "
                          "convergence may be slower");
  if (c.firstIterationRescale || c.asymmetricErrors)
    hf_errlog(2026102213, "W: GaussNewton: FirstIterationRescale or asymmetric errors make the "
                          "derivative model approximate, convergence may be slower");
}

DataSnapshot FortranBridge::snapshot() const {
  const int n = _config.fitPoints, m = _config.nuisances;
  DataSnapshot s;
  s.data.resize(n);
  s.theory.resize(n);
  s.sigma.resize(n);
  s.sigma2DT.resize(n);
  s.sigma2TT.resize(n);
  s.alpha.resize(n);
  s.gamma.resize(m, n);  // column-major, as the Fortran Gamma(NNuis, NFit)
  s.gammaPower.resize(m);
  s.shifts.resize(m);
  s.priorScale.resize(m);
  std::vector<int> additive(n), logTerm(n);
  gn_bridge_snapshot_(n, _points.data(), m, _sources.data(), s.data.data(), s.theory.data(),
                      s.sigma.data(), s.sigma2DT.data(), s.sigma2TT.data(), additive.data(),
                      logTerm.data(), s.alpha.data(), s.gamma.data(), s.gammaPower.data(),
                      s.shifts.data(), s.priorScale.data());
  s.additive.assign(additive.begin(), additive.end());
  s.logTerm.assign(logTerm.begin(), logTerm.end());
  if (_config.extraSystRescale) {
    s.extraSystRescale = true;
    s.sigmaFirstPass.resize(n);
    s.sigma2DTStat.resize(n);
    s.betaLinear.resize(m, n);
    gn_bridge_rescale_(n, _points.data(), m, _sources.data(), s.sigmaFirstPass.data(), s.sigma2DTStat.data(),
                       s.betaLinear.data());
  }
  return s;
}

Eigen::VectorXd FortranBridge::theory() const {
  Eigen::VectorXd t(_config.fitPoints);
  gn_bridge_theory_(_config.fitPoints, _points.data(), t.data());
  return t;
}

double FortranBridge::polarisationChi2() const {
  double chi2 = 0;
  gn_bridge_polarisation_chi2_(chi2);
  return chi2;
}

}
}
