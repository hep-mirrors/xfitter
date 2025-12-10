#pragma once

#include "BaseEvolution.h"
#include"BasePdfDecomposition.h"

#include <apfel/grid.h>
#include <apfel/dglap.h>
#include <apfel/dglapbuilder.h>
#include <apfel/tabulateobject.h>
#include <apfel/alphaqcdmsbarmass.h>

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
    ///@}

  private:
    BasePdfDecomposition*                                                   _inPDFs;
    std::vector<double>                                                     _Masses;
    std::vector<double>                                                     _Thresholds;
    std::unique_ptr<const apfel::Grid>                                      _Grid;
    std::map<int,apfel::DglapObjects>                                       _DglapObj;
    std::unique_ptr<apfel::TabulateObject<apfel::Set<apfel::Distribution>>> _TabulatedPDFs;
    std::function<double(double const& Q)>                                  _AlphaQCD;
    /// pertirbative order
    int _PtOrder;
    /// resummation scale coefficient
    double _xi;
    /// ratio of resummation scale in alphaS evolution over resummation scale in PDF evolution
    double _xi_aspdf;
    /// pointer to alphas parameter
    double* _alphas;
    /// Evolution starting scale:
    double _Q0;
    /// Evolution starting scale for alphas:
    const double* _alphas_q0;
    /// Flavour scheme
    int _isFFNS;
    /// Number of flavours
    int _NFlavour;
    /// heavy quark mass scheme (pole or MSbar)
    std::string _heavyQuarkMassScheme;
    /// Heavy-quark masses
    const double *_mch;
    const double *_mbt;
    const double *_mtp;
    /// Heavy-quark masses at last iteration
    double _mch_last;
    double _mbt_last;
    double _mtp_last;
  };
}
