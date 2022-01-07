
#pragma once

#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"
#include <memory>

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
    UvDvUbarDbarS(const char*name);
    virtual const char*getClassName()const override final;

    /// Optional initialization at the first call
    virtual void atStart()override final;

    /// Compute sum-rules
    virtual void atIteration() override final;
    
    /// Compute PDF in a physical base in LHAPDF format at the initial scale
    virtual std::map<int,double>xfxMap(double x)const override final;
  private:
    BasePdfParam*par_xuv{nullptr},
      *par_xdv{nullptr},
      *par_xubar{nullptr},
      *par_xdbar{nullptr},
      *par_xs{nullptr},
      *par_xg{nullptr},
      *par_xgamma{nullptr};
  };
}
