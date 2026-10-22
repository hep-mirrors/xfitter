#include "DerivativeCheck.h"
#include "ProfiledChi2.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace xfitter {
namespace gn {

double checkDerivatives(ProfiledChi2& problem, const Eigen::VectorXd& p, const Eigen::VectorXd& steps,
                        const std::vector<std::string>& names) {
  const double chi2 = problem.chi2(p);
  const QuadraticModel model = problem.model(p);

  std::printf("\nGaussNewton derivative check at the start point, chi2 = %.6f\n", chi2);
  std::printf("%5s %15s %15s %15s %12s\n", "", "parameter", "model dchi2/dp", "FD dchi2/dp", "rel. diff");
  double worst = 0;
  for (int j = 0; j < p.size(); ++j) {
    Eigen::VectorXd up = p, down = p;
    up(j) += steps(j);
    down(j) -= steps(j);
    const double fd = (problem.chi2(up) - problem.chi2(down)) / (2 * steps(j));
    const double rel = std::abs(model.gradient(j) - fd) / std::max(std::abs(fd), 1e-12);
    worst = std::max(worst, rel);
    std::printf("%5d %15s %15.6g %15.6g %12.2e\n", j, names[j].c_str(), model.gradient(j), fd, rel);
  }

  // Hessian along a few directions: second difference of the official chi2
  std::srand(20261022);
  for (int k = 0; k < 3; ++k) {
    const Eigen::VectorXd v = Eigen::VectorXd::Random(p.size()).cwiseProduct(steps);
    const double fd = problem.chi2(p + v) - 2 * chi2 + problem.chi2(p - v);
    const double gn = v.dot(model.hessian * v);
    std::printf("direction %d: v.H.v model = %.6g, second difference = %.6g, rel. diff = %.2e\n", k, gn, fd,
                std::abs(gn - fd) / std::abs(fd));
  }
  std::printf("GaussNewton derivative check: largest relative gradient deviation %.2e\n\n", worst);
  return worst;
}

}
}
