#include "bridge/InputPdfs.h"
#include "bridge/TheoryEvaluator.h"
#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"
#include "dependent_pars.h"
#include "xfitter_pars.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace xfitter {
namespace gn {

namespace {
const int kFlavours[] = {-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6, 22};  // 22: photon (QED)
const int kNumFlavours = sizeof(kFlavours) / sizeof(kFlavours[0]);

/// 100 nodes uniform in log x on [1e-6, 0.1), 50 uniform in x on [0.1, 0.98]
std::vector<double> nodes() {
  std::vector<double> x;
  for (int i = 0; i < 100; ++i) x.push_back(std::pow(10., -6. + 5. * i / 100));
  for (int i = 0; i < 50; ++i) x.push_back(0.1 + 0.88 * i / 49);
  return x;
}
}  // namespace

InputPdfs::InputPdfs(TheoryEvaluator& evaluator) : _evaluator(evaluator), _x(nodes()) {}

int InputPdfs::size() const {
  return static_cast<int>(XFITTER_PARS::gPdfDecompositions.size() * kNumFlavours * _x.size());
}

Eigen::VectorXd InputPdfs::sample(const Eigen::VectorXd& p) {
  // the first part of init_at_iteration (src/ftheor_eval.cc)
  _evaluator.setParameters(p);
  updateDependentParameters();
  for (const auto& it : XFITTER_PARS::gParameterisations) it.second->atIteration();
  for (const auto& it : XFITTER_PARS::gPdfDecompositions) it.second->atIteration();

  Eigen::VectorXd f(size());
  int k = 0;
  for (const auto& it : XFITTER_PARS::gPdfDecompositions)
    for (const double x : _x) {
      const std::map<int, double> xf = it.second->xfxMap(x);
      for (const int flavour : kFlavours) {
        const auto found = xf.find(flavour);
        f(k++) = found != xf.end() ? found->second : 0.;
      }
    }
  if (_weights.size() == 0 && f.allFinite()) setWeights(f);
  return _weights.size() ? Eigen::VectorXd(f.cwiseProduct(_weights)) : f;
}

void InputPdfs::setWeights(const Eigen::VectorXd& f) {
  // each node counts by its relative change, as the data do; plain xf values are dominated
  // by the smallest x for a steep gluon (D13)
  const int nx = static_cast<int>(_x.size());
  const double global = f.cwiseAbs().maxCoeff();
  _weights.resize(f.size());
  for (int offset = 0; offset < f.size(); offset += kNumFlavours * nx)
    for (int flavour = 0; flavour < kNumFlavours; ++flavour) {
      double largest = 0;
      for (int a = 0; a < nx; ++a) largest = std::max(largest, std::abs(f(offset + a * kNumFlavours + flavour)));
      for (int a = 0; a < nx; ++a) {
        const int k = offset + a * kNumFlavours + flavour;
        // a flavour that vanishes at the first point (c, b, t at Q0) is measured in global units
        const double scale = largest > 0 ? std::max(std::abs(f(k)), 1e-3 * largest) : global;
        _weights(k) = scale > 0 ? 1 / scale : 1.;
      }
    }
}

}
}
