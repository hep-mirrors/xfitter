#pragma once

// J_t by forward or central finite differences of the theory, without chi2 or nuisance
// solve per column (README.md, Sec. 4, "Design decisions: FD step").

#include "jacobian/TheoryJacobian.h"

#include <vector>

namespace xfitter {
namespace gn {
class TheoryEvaluator;

struct FDOptions {
  bool central = false;
  double relativeStep = 1e-3;  ///< h_j = relativeStep * max(|p_j|, yaml step_j)
  int workers = 0;             ///< forked processes for the columns, 0: none
};

class FiniteDifferenceJacobian : public TheoryJacobian {
public:
  FiniteDifferenceJacobian(TheoryEvaluator& evaluator, const Eigen::VectorXd& yamlSteps,
                           const Eigen::VectorXd& lower, const Eigen::VectorXd& upper,
                           const FDOptions& options);

  Eigen::MatrixXd compute(const Eigen::VectorXd& p, const Eigen::VectorXd& t0) override;
  int evaluations() const override { return _evaluations; }

  /// The columns `which` of J_t only (n points x which.size())
  Eigen::MatrixXd columns(const Eigen::VectorXd& p, const Eigen::VectorXd& t0, const std::vector<int>& which);

  /// Column j is (t(p + up e_j) - t(p + down e_j)) / (up - down); {0, 0} if there is no room
  struct Displacement {
    double up = 0, down = 0;
  };
  Displacement displacement(const Eigen::VectorXd& p, int j) const {
    return displacement(p, j, _options.central);
  }
  /// The same step h_j, central if `central` and both sides have room
  Displacement displacement(const Eigen::VectorXd& p, int j, bool central) const;

private:
  double step(const Eigen::VectorXd& p, int j) const;
  /// Column j: one or two theory evaluations, respecting the bounds
  Eigen::VectorXd column(const Eigen::VectorXd& p, const Eigen::VectorXd& t0, int j);

  TheoryEvaluator& _evaluator;
  Eigen::VectorXd _yamlSteps, _lower, _upper;
  FDOptions _options;
  int _evaluations = 0;
};

}
}
