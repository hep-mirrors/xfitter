/*
   @file GaussNewtonMinimizer.cc
   @date 2026-10-22

   xFitter glue of the GaussNewton minimiser. See README.md.
*/

#include "GaussNewtonMinimizer.h"
#include "xfitter_cpp_base.h"

#include <cmath>

namespace xfitter {

/// the class factory, for dynamic loading
extern "C" GaussNewtonMinimizer* create() {
  return new GaussNewtonMinimizer();
}

GaussNewtonMinimizer::GaussNewtonMinimizer() : BaseMinimizer("GaussNewton") {
}

void GaussNewtonMinimizer::atStart() {
}

void GaussNewtonMinimizer::doMinimization() {
  hf_errlog(2026102210, "F: GaussNewton minimizer is not implemented yet");
}

void GaussNewtonMinimizer::errorAnalysis() {
}

ConvergenceStatus GaussNewtonMinimizer::convergenceStatus() {
  return _convergenceStatus;
}

double GaussNewtonMinimizer::getParameterUncertainty(const std::string& name) const {
  int index = parameterIndex(name);
  if (index < 0 || static_cast<size_t>(index) >= _parameterUncertainties.size()) return std::nan("");
  return _parameterUncertainties[index];
}

}
