#include "FKTables.h"

#include <algorithm>

namespace xfitter {
namespace fk {

OperatorTable TabulateFK(StructureFunctionObjectsFunc const& FObj, std::vector<int> const& components,
                         int PerturbativeOrder, std::function<double(double const&)> const& Alphas,
                         std::function<std::vector<double>(double const&)> const& Couplings,
                         OperatorTable const& evolution, int nQ, double QMin, double QMax, int InterDegree,
                         std::vector<double> const& Thresholds) {
  const auto fkAtQ = [=, &evolution] (double const& Q) -> apfel::Set<apfel::Operator> {
    const apfel::StructureFunctionObjects FObjQ = FObj(Q, Couplings(Q));
    const apfel::Set<apfel::Operator> E = evolution.Evaluate(Q);
    const apfel::Operator zero = 0. * E.at(0);
    // coefficient functions summed over the perturbative orders, as BuildStructureFunctions
    // (xiR = xiF = 1)
    const double cp = Alphas(Q) / apfel::FourPi;
    std::map<int, apfel::Operator> fk;
    for (int j = 0; j <= 12; ++j) fk.insert({j, zero});
    for (const int k : components) {
      apfel::Set<apfel::Operator> C = FObjQ.C0.at(k);
      if (PerturbativeOrder > 0) C += cp * FObjQ.C1.at(k);
      if (PerturbativeOrder > 1) C += cp * cp * FObjQ.C2.at(k);
      if (PerturbativeOrder > 2) C += cp * cp * cp * FObjQ.C3.at(k);
      for (int j = 0; j <= 12; ++j) {
        // column j of the evolution operator in place of the evolved PDFs; the distributions
        // in FObjQ.skip are set to zero, as in BuildStructureFunctions
        std::map<int, apfel::Operator> column;
        bool any = false;
        for (int i = 0; i <= 12; ++i) {
          const auto it = apfel::Gkj.find({i, j});
          const bool skipped = std::find(FObjQ.skip.begin(), FObjQ.skip.end(), i) != FObjQ.skip.end();
          if (it != apfel::Gkj.end() && !skipped) {
            column.insert({i, E.at(it->second)});
            any = true;
          } else {
            column.insert({i, zero});
          }
        }
        if (any) fk.at(j) += (C * apfel::Set<apfel::Operator>{FObjQ.ConvBasis.at(k), column}).Combine();
      }
    }
    return apfel::Set<apfel::Operator>{apfel::EvolveDistributionsBasisQCD{}, fk};
  };
  return OperatorTable{fkAtQ, nQ, QMin, QMax, InterDegree, Thresholds};
}

void Row::add(OperatorTable const& table, double x, double Q, double factor) {
  if (factor == 0) return;
  const apfel::Set<apfel::Operator> fk = table.Evaluate(Q);
  for (int j = 0; j <= 12; ++j) {
    const std::vector<double>& w = fk.at(j).Evaluate(x).GetDistributionJointGrid();
    double* out = _w.data() + static_cast<size_t>(j) * _nJoint;
    for (int a = 0; a < _nJoint; ++a) out[a] += factor * w[a];
  }
}

double Row::contract(std::vector<double> const& f0) const {
  double s = 0;
  for (size_t i = 0; i < _w.size(); ++i) s += _w[i] * f0[i];
  return s;
}

std::vector<double> Flatten(std::map<int, apfel::Distribution> const& f0) {
  const int nJoint = static_cast<int>(f0.at(0).GetDistributionJointGrid().size());
  std::vector<double> out(13 * static_cast<size_t>(nJoint), 0.);
  for (auto const& d : f0) {
    if (d.first < 0 || d.first > 12) continue;
    const std::vector<double>& v = d.second.GetDistributionJointGrid();
    std::copy(v.begin(), v.end(), out.begin() + static_cast<size_t>(d.first) * nJoint);
  }
  return out;
}

}
}
