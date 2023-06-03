#pragma once

/**
  @class' EvolutionLHAPDF

  @brief A wrapper class for LHAPDF evolution

  @version 0.1
  @date 2018-08-20
  */

#include "BaseEvolution.h"
//Try to suppress unused-local-typedef warning from boost 1.53.0 for gcc
//Apparently these warnings have been fixed in later versions of boost
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "LHAPDF/LHAPDF.h"
#pragma GCC diagnostic pop

namespace xfitter {

  class EvolutionLHAPDF : BaseEvolution {
  public:
    /// Empty constructor (needed for the dynamic loading)
    EvolutionLHAPDF(const char* name);
    virtual const char* getClassName()const override final;

  public:
    /// Global initialization
    virtual void atStart() override final;
    /// Init at each change of at least one parameter
    virtual void atConfigurationChange() override final;
    virtual std::map<int, double>xfxQmap(double x, double Q)override final;
    virtual double xfxQ(int i, double x, double Q)override final;
    virtual void xfxQarray(double x, double Q, double* pdfs)override final;
    virtual double getAlphaS(double Q)override final;
    /// Get property
    virtual std::string getPropertyS(std::string const& propertyName) const override final;

    /// Get property
    virtual int  getPropertyI(std::string const& propertyName) const override final;

  /// Get property
    virtual double getPropertyD(std::string const& propertyName) const override final;

  /// Get property with default value
    virtual double getPropertyD(std::string const& propertyName, double defval) const override final;
    
  private:
    std::string _set_name{""};
    int _member{0};
    LHAPDF::PDF* _pdf{nullptr};
    bool _has_photon;
  };
};  // namespace xfitter
