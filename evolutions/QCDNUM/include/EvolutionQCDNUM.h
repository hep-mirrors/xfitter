#pragma once

#include "BaseEvolution.h"
#include <vector>
#include <memory>

/**
  @class' EvolutionQCDNUM

  @brief A wrapper class for QCDNUM evolution

  @version 0.2
  @date 2018-09-29
  */

namespace xfitter
{

  class EvolutionQCDNUM: public BaseEvolution
  {
  public:
    EvolutionQCDNUM(const char* name):BaseEvolution{name}{};

    virtual const char* getClassName() const override final;
    virtual void atStart() override final;
    virtual void atIteration() override final;
    virtual void afterIteration() override final;
    virtual void atConfigurationChange() override final;
    virtual std::map<int,double>xfxQmap(double x, double Q) override final;
    virtual double xfxQ(int i, double x, double Q) override final;
    virtual void xfxQarray(double x, double Q, double*pdfs) override final;
    virtual double getAlphaS(double Q) override final;
    virtual std::vector<double> getXgrid() override final;
    virtual std::vector<double> getQgrid() override final;

    /// Helper to get PDF type
    const int getPdfType() const {return _itype;}

  private:
    /// PDFs called outside boundaries check:
    int _icheck{0};

    /// Number of additional to quark PDFs (e.g. photon)
    int _nExt{0};

    /// PDF type (1 -- unpolorized internal)
    int _itype{1};

    /// Spline order for interpolation
    int _splineOrder{2};

    /// Read or not QCDNUM pre-stored tables (makes init faster)
    int _readTables{0};

    const double* Mz;/// Z-boson mass
    const double* alphas;/// alphaS(Mz^2)
  };
};
