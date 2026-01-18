#pragma once

#include "BaseEvolution.h"
#include"BasePdfDecomposition.h"

#include <vector>
#include <memory>

namespace xfitter
{
  /**
     @class EvolutionOPENQCDRAD

     @brief Derived class of BaseEvolution for using OPENQCDRAD as an evolution code.

     @version 0.1
     @date 2026-01-18
  */
  class EvolutionOPENQCDRAD: public BaseEvolution
  {
  public:
    EvolutionOPENQCDRAD(const char*name):BaseEvolution{name}{}
    virtual const char*getClassName()const final override;
    virtual void atStart()override final;
    virtual void atIteration()override final;

    /**
     * @name Getters
     */
    virtual std::map<int,double>xfxQmap(double x,double Q)override final;
    virtual double xfxQ(int i,double x,double Q)override final;
    virtual void xfxQarray(double x,double Q,double*pdfs)override final;
    virtual double getAlphaS(double Q)override final;
    //virtual std::vector<double> getXgrid() override final;
    //virtual std::vector<double> getQgrid() override final;
    ///@}

  private:
    int _order;
    //double _Qmin;
    //double _Qmax;
    //double _ymax;
    //int _isFFNS;
    int _nflavour;
    int _msbar;
    const double* _mch;
    const double* _mbt;
    //const double* _mtp;
    /// pointer to alphas parameter
    double* _alphas;
    /// Evolution starting scale for alphas:
    const double* _alphas_q0;
    BasePdfDecomposition* _inPDFs;
  };
}
