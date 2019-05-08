#pragma once
#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"

/**
  @class UvDvUbarDbarSSbarPdfDecomposition

  @brief A class for UvDvUbarDbarSSbar pdf decomposition

  @version 0.1
  @date 2019-02-25
  */

namespace xfitter {

  class UvDvUbarDbarSSbarPdfDecomposition : public BasePdfDecomposition
  {
  public:
    /// Default constructor. Name is the PDF name
    UvDvUbarDbarSSbarPdfDecomposition (const char *inName);

    virtual const char*getClassName()const override final;

    /// Optional initialization at the first call
    virtual void atStart() override final;

    /// Compute sum-rules
    virtual void atIteration() override final;

    /// Compute PDF in a physical base in LHAPDF format for given x and Q
    virtual std::map<int,double>xfxMap(double x)const override final;

  private:
    BasePdfParam*par_xuv{nullptr},
    *par_xdv{nullptr},
    *par_xubar{nullptr},
    *par_xdbar{nullptr},
    *par_xs{nullptr},
    *par_xsbar{nullptr},
    *par_xg{nullptr};
  };
}
