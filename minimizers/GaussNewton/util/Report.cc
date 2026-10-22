#include "util/Report.h"
#include "xfitter_cpp.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace xfitter {
namespace gn {

namespace {

std::string outputFile(const std::string& name) {
  return stringFromFortran(coutdirname_.outdirname, sizeof(coutdirname_.outdirname)) + "/" + name;
}

std::string fixed6(double x) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%.6f", x);
  return buffer;
}

void writeMatrix(std::ostream& out, const std::vector<std::string>& names, const Eigen::MatrixXd& M) {
  out << std::setw(14) << " ";
  for (const std::string& name : names) out << std::setw(14) << name;
  out << std::endl << std::setprecision(4);
  for (int i = 0; i < M.rows(); ++i) {
    out << std::setw(14) << names[i];
    for (int j = 0; j < M.cols(); ++j) out << std::setw(14) << M(i, j);
    out << std::endl;
  }
}

}  // namespace

void writeParameters(const std::vector<std::string>& names, const Eigen::VectorXd& values,
                     const Eigen::VectorXd& errors) {
  std::ofstream parsout(outputFile("parsout_0"));
  if (!parsout) hf_errlog(2026102240, "W: GaussNewton: failed to open parsout_0");
  for (size_t i = 0; i < names.size(); ++i)
    parsout << std::setw(5) << i << "   '" << names[i] << "'    " << values(i) << "    " << errors(i) << std::endl;

  std::ofstream yaml(outputFile("pars.yaml"));
  yaml << "Parameters:" << std::endl;
  for (size_t i = 0; i < names.size(); ++i)
    yaml << "  " << names[i] << ": [ " << fixed6(values(i)) << ", " << fixed6(errors(i)) << " ]" << std::endl;
}

void writeCovariance(const std::vector<std::string>& names, const Eigen::MatrixXd& covariance,
                     const std::vector<std::string>& summary) {
  std::ofstream out(outputFile("gaussnewton.out.txt"));
  if (!out) hf_errlog(2026102241, "W: GaussNewton: failed to open gaussnewton.out.txt");
  for (const std::string& line : summary) out << line << std::endl;

  const Eigen::VectorXd sigma = covariance.diagonal().cwiseSqrt();
  const Eigen::MatrixXd correlation = sigma.cwiseInverse().asDiagonal() * covariance * sigma.cwiseInverse().asDiagonal();
  out << std::endl << "COVARIANCE MATRIX " << std::endl;
  writeMatrix(out, names, covariance);
  out << std::endl << "CORRELATION MATRIX " << std::endl;
  writeMatrix(out, names, correlation);

  const Eigen::MatrixXd inverse = covariance.inverse();
  out << "GLOBAL CORRELATION COEFFICIENT" << std::endl;
  for (size_t i = 0; i < names.size(); ++i)
    out << std::setw(5) << i << std::setw(15) << names[i] << std::setw(15)
        << std::sqrt(1. - 1. / (covariance(i, i) * inverse(i, i))) << std::endl;
}

void printParameters(const std::vector<std::string>& names, const Eigen::VectorXd& values,
                     const Eigen::VectorXd& errors) {
  std::cout << std::endl << "Fitted parameters" << std::endl;
  for (size_t i = 0; i < names.size(); ++i)
    std::cout << std::setw(5) << i << std::setw(15) << names[i] << std::setw(15) << values(i)
              << " +/- " << errors(i) << std::endl;
  std::cout << std::endl << "----- Parameters in YAML format (can copy-paste into parameters.yaml):" << std::endl;
  std::cout << "  Parameters:" << std::endl;
  for (size_t i = 0; i < names.size(); ++i)
    std::cout << "  " << names[i] << ": [ " << fixed6(values(i)) << ", " << fixed6(errors(i)) << " ]" << std::endl;
  std::cout << "----- End of parameters in YAML format" << std::endl << std::endl;
}

}
}
