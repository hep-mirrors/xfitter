#pragma once

// Output files of the GaussNewton minimizer, in the formats used by CERES:
// parsout_0, pars.yaml and gaussnewton.out.txt (covariance, correlation, summary).

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace xfitter {
namespace gn {

/// parsout_0 and pars.yaml in the xFitter output directory
void writeParameters(const std::vector<std::string>& names, const Eigen::VectorXd& values,
                     const Eigen::VectorXd& errors);

/// gaussnewton.out.txt: the summary lines, then covariance, correlation and global
/// correlation coefficients
void writeCovariance(const std::vector<std::string>& names, const Eigen::MatrixXd& covariance,
                     const std::vector<std::string>& summary);

/// Parameters, errors and YAML block on stdout
void printParameters(const std::vector<std::string>& names, const Eigen::VectorXd& values,
                     const Eigen::VectorXd& errors);

}
}
