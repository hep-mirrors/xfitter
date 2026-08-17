
#pragma once

#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"
#include <memory>
#include <string>

namespace xfitter
{
  class MSHT_Decomposition : public BasePdfDecomposition
  {
  public:
    MSHT_Decomposition(const char* name);
    virtual const char* getClassName() const override final;

    virtual void atStart() override final;
    virtual void atIteration() override final;
    virtual std::map<int,double> xfxMap(double x) const override final;

  private:
    BasePdfParam *par_xuv{nullptr},
                 *par_xdv{nullptr},
                 *par_S{nullptr},
                 *par_splus{nullptr},
                 *par_sminus{nullptr},
                 *par_d_over_u{nullptr},
                 *par_xg{nullptr};
  };
}