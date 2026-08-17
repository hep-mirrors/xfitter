#pragma once
#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"

/**
  @class MSHT20PdfDecomposition

  @brief A class for MSHT20 pdf decomposition

  @version 0.1
  @date 2026-08-16
  */

namespace xfitter {

  class MSHT20PdfDecomposition : public BasePdfDecomposition
  {
  public:
    /// Default constructor. Name is the PDF name
    MSHT20PdfDecomposition (const char *inName);

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
    *par_xs{nullptr},
    *par_xr{nullptr},
    *par_xsp{nullptr},
    *par_xsm{nullptr},
    *par_xg{nullptr},
    *par_xgn{nullptr};
  };
}
