#pragma once

#include "BaseEvolution.h"
#include"BasePdfDecomposition.h"

#include <vector>
#include <memory>
#include <functional>

namespace xfitter
{
  /**
     @class EvolutionAPFELxx

     @brief Derived class of BaseEvolution for using APFEL++ as an evolution code.

     @version 0.2
     @date 2018-09-29
  */
  class EvolutionHoppet: public BaseEvolution
  {
  public:
    EvolutionHoppet(const char*name):BaseEvolution{name}{}
    virtual const char*getClassName()const final override;

    /**
     * @brief Function that initialises the evolution in APFEL++.
     */
    virtual void atStart()override final;
    virtual void atIteration()override final;

    /**
     * @name Getters
     */
    virtual std::map<int,double>xfxQmap(double x,double Q)override final;
    virtual double xfxQ(int i,double x,double Q)override final;
    virtual void xfxQarray(double x,double Q,double*pdfs)override final;
    virtual double getAlphaS(double Q)override final;
    virtual std::vector<double> getXgrid() override final;
    virtual std::vector<double> getQgrid() override final;
    ///@}

  private:
    BasePdfDecomposition*                                                   _inPDFs;
    //std::vector<double>                                                     _Masses;
    //std::vector<double>                                                     _Thresholds;
    //std::function<double(double const& Q)>                                  _AlphaQCD;
  };
  
  

}
