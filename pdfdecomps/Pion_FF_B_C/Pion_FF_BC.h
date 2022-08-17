
#pragma once

#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"
#include <memory>

namespace xfitter
{
  /**
     @class Pion_FF_BC_PdfDecomposition 

     @brief A class for Pion-FF_ pdf decomposition

     @version 0.1
     @date 2018-07-11
  */

  class Pion_FF_BC : public BasePdfDecomposition
  {
  public:
    /// Default constructor. Name is the PDF name
    Pion_FF_BC(const char*name);
    virtual const char*getClassName()const override final;

    /// Optional initialization at the first call
    virtual void atStart()override final;

    /// Compute sum-rules
    virtual void atIteration() override final;
    
    /// Compute PDF in a physical base in LHAPDF format at the initial scale
    virtual std::map<int,double>xfxMap(double x)const override final;
  private:
    BasePdfParam*par_xup{nullptr},
                *par_xcp{nullptr},
                *par_xbp{nullptr},
                *par_xsp{nullptr},
                *par_xg{nullptr};
  };
}
