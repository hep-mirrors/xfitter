
#pragma once

#include "BasePdfDecomposition.h"
//Try to suppress unused-local-typedef warning from boost 1.53.0 for gcc
//Apparently these warnings have been fixed in later versions of boost
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "LHAPDF/LHAPDF.h"
#pragma GCC diagnostic pop


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
    LHAPDFDecomposition(const char*name);
    ~LHAPDFDecomposition();
    virtual const char*getClassName()const override final;

    /// Optional initialization at the first call
    virtual void atStart()override final;

    /// Compute PDF in a physical basis in LHAPDF format at the initial scale
    virtual std::map<int,double>xfxMap(double x)const override final;

  private:
    LHAPDF::PDF*_pdf{nullptr};
  };
}
