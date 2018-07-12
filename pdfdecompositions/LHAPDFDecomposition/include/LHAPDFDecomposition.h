
#pragma once

#include "BasePdfDecomposition.h"

#include "LHAPDF/LHAPDF.h"

/**
  @class LHAPDFDecomposition

  @brief A class for LHAPDF pdf decomposition

  @version 0.1
  @date 2018-07-12
  */

namespace xfitter
{
  class LHAPDFDecomposition : public BasePdfDecomposition
  {
  public:
    /// Default constructor. Name is the PDF name
    LHAPDFDecomposition(const std::string& PDFset, const int& mem = 0);

    /// Optional initialization at the first call
    virtual void initAtStart(const std::string& pars) const override final;

    /// Compute PDF in a physical base in LHAPDF format at the initial scale
    virtual std::function<std::map<int,double>(const double& x)> f0() const override final;

  private:
    const int                 _mem;
    std::vector<LHAPDF::PDF*> _dist;
  };
}
