#pragma once

#include "BaseEvolution.h"
#include"BasePdfDecomposition.h"

#include <apfel/grid.h>
#include <apfel/dglap.h>
#include <apfel/dglapbuilder.h>
#include <apfel/tabulateobject.h>

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
  class EvolutionAPFELxx: public BaseEvolution
  {
  public:
    EvolutionAPFELxx(const char*name):BaseEvolution{name}{}
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

    apfel::TabulateObject<apfel::Set<apfel::Distribution>> GetTabulatedPDFs();
    std::function<double(double const& Q)>                 GetAlphaQCD();

    /**
     * @brief PDFs and alphas in a fixed-flavour-number scheme with
     * nf = 3 or 4 active flavours (nf = 5 returns the nominal VFNS
     * PDFs and alphas)
     */
    std::map<int,double> xfxQmapNF(int nf,double x,double Q);
    double getAlphaSNF(int nf,double Q);
    ///@}

  private:
    BasePdfDecomposition*                                                   _inPDFs;
    std::vector<double>                                                     _Masses;
    std::vector<double>                                                     _Thresholds;
    std::unique_ptr<const apfel::Grid>                                      _Grid;
    std::map<int,apfel::DglapObjects>                                       _DglapObj;
    std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> _TabulatedPDFs;
    std::function<double(double const& Q)>                                  _AlphaQCD;

    // 3- and 4-flavour schemes
    std::vector<double>                                                     _Masses3F;
    std::vector<double>                                                     _Thresholds3F;
    std::map<int,apfel::DglapObjects>                                       _DglapObj3F;
    std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> _TabulatedPDFs3F;
    std::function<double(double const& Q)>                                  _AlphaQCD3F;

    std::vector<double>                                                     _Masses4F;
    std::vector<double>                                                     _Thresholds4F;
    std::map<int,apfel::DglapObjects>                                       _DglapObj4F;
    std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> _TabulatedPDFs4F;
    std::function<double(double const& Q)>                                  _AlphaQCD4F;
  };
}
