
#pragma once

#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"

namespace xfitter
{
  /**
     @class UvDvUbarDbarS_PdfDecomposition 

     @brief A class for UvDvUbarDbarS_ pdf decomposition

     @version 0.1
     @date 2018-07-11
  */

  class UvDvUbarDbarS : public BasePdfDecomposition
  {
  public:
    /// Default constructor. Name is the PDF name
    UvDvUbarDbarS();

    /// Optional initialization at the first call
    virtual void initAtStart(const std::string& pars) override final;

    /// Compute PDF in a physical base in LHAPDF format at the initial scale
    virtual std::function<std::map<int,double>(const double& x)> f0() const override final;

  private:
    std::map<std::string,BasePdfParam*> pdf_pars;
  };
}
