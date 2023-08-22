
#pragma once

#include "BasePdfDecomposition.h"
#include "BaseEvolution.h"


/**
  @class FromEvolutionDecomposition

  @brief A class for decomposition which is based on another PDF evolution (need to supply the evolution name and the starting scale)

  @version 0.1
  @date 2023-08-18
  */

namespace xfitter
{
  class FromEvolutionDecomposition : public BasePdfDecomposition
  {
  public:
    FromEvolutionDecomposition(const char*name);
    ~FromEvolutionDecomposition();
    virtual const char*getClassName()const override final;

    /// Optional initialization at the first call
    virtual void atStart()override final;

    /// Optional initialization at each iteration. Can be used to compute sum-rules
    virtual void atIteration()override final;

    /// Compute PDF in a physical basis in LHAPDF format at the initial scale
    virtual std::map<int,double>xfxMap(double x)const override final;

  private:
    BaseEvolution* _evolution{nullptr};
    //std::string _Q0_parname;
    double _Q0;
  };
}
