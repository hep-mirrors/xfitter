#pragma once

// The input PDFs f0(p): xf(x) at Q0 of every decomposition, sampled at fixed x nodes and
// measured relative to the first sample (README.md, Eq. 14 and D13). Only the
// parameterisation part of an iteration runs
// (dependent parameters, parameterisations, decompositions with the sum rules): no
// evolution, no theory. Changes the global xFitter state like a theory evaluation.

#include <Eigen/Dense>

#include <vector>

namespace xfitter {
namespace gn {
class TheoryEvaluator;

class InputPdfs {
public:
  explicit InputPdfs(TheoryEvaluator& evaluator);

  /// f0(p): for each decomposition and node, xf of the flavours -6..6 and the photon, times
  /// 1 / max(|xf|, 1e-3 max_x |xf|) of the first finite sample (fixed afterwards)
  Eigen::VectorXd sample(const Eigen::VectorXd& p);

  int size() const;

private:
  TheoryEvaluator& _evaluator;
  void setWeights(const Eigen::VectorXd& f);

  std::vector<double> _x;    ///< nodes
  Eigen::VectorXd _weights;  ///< empty until the first finite sample
};

}
}
