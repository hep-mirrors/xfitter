#pragma once

// FK tables of DIS structure functions from APFEL++ operators.
//
// At fixed alpha_s, masses and Q0 a structure function is linear in the input PDFs at Q0:
//   F(x, Q) = sum_j sum_a K_{j,a}(x, Q) f0_j(x_a)
// with j the QCD evolution basis index and x_a the nodes of the evolution grid. The FK
// operator of input j is C(Q) * E_{.j}(Q): the coefficient functions times the column j of
// the DGLAP evolution operator, built exactly as BuildStructureFunctions builds the
// structure functions from the PDFs (same convolution basis, same "skip" list), and
// tabulated in Q like the structure functions of the reactions. The row K(x, Q) is the
// interpolated operator evaluated at x.

#include <apfel/apfelxx.h>

#include <functional>
#include <map>
#include <vector>

namespace xfitter {
namespace fk {

using StructureFunctionObjectsFunc =
    std::function<apfel::StructureFunctionObjects(double const&, std::vector<double> const&)>;
using OperatorTable = apfel::TabulateObject<apfel::Set<apfel::Operator>>;

/// Tabulated FK operators of one structure-function piece: Q -> {input j -> operator}.
/// components: structure-function indices k summed (e.g. {1,2,3} for the light part);
/// evolution: tabulated DGLAP evolution operator of the flavour scheme of FObj.
OperatorTable TabulateFK(StructureFunctionObjectsFunc const& FObj, std::vector<int> const& components,
                         int PerturbativeOrder, std::function<double(double const&)> const& Alphas,
                         std::function<std::vector<double>(double const&)> const& Couplings,
                         OperatorTable const& evolution, int nQ, double QMin, double QMax, int InterDegree,
                         std::vector<double> const& Thresholds);

/// FK row: weights per input j (0..12) and joint-grid node, flattened as j * nJoint + a
class Row {
public:
  explicit Row(int nJoint = 0) : _nJoint(nJoint), _w(13 * static_cast<size_t>(nJoint), 0.) {}
  /// Adds factor * (row of table at Q, evaluated at x)
  void add(OperatorTable const& table, double x, double Q, double factor);
  /// sum_j sum_a w_{j,a} f0_j(x_a)
  double contract(std::vector<double> const& f0) const;

private:
  int _nJoint;
  std::vector<double> _w;
};

/// Input PDFs at Q0 on the joint grid, flattened as j * nJoint + a (QCD evolution basis)
std::vector<double> Flatten(std::map<int, apfel::Distribution> const& f0);

}
}
