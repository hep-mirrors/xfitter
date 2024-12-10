#pragma once

#include "BaseEvolution.h"
#include <vector>
#include <memory>

/**
  @class' EvolutionPBTMD

  @brief A wrapper class for PBTMD evolution

  @version 0.2
  @date 2018-09-29
  */

namespace xfitter
{

  class EvolutionPBTMD: public BaseEvolution
  {
  public:
//    EvolutionPBTMD(const char* name):BaseEvolution{name}{};
    /// Empty constructor (needed for the dynamic loading)
    EvolutionPBTMD(const char* name);
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
    virtual void Write_TMD(const char* name) override final;

    /// Helper to get PDF type
    const int getPdfType() const {return _itype;}

  private:
    std::string _set_name{""};
    int _member{0};
    /// PDF type (1 -- unpolorized internal)
    int _itype{1};
    double xSav{-1}, Q2Sav{-1.} ;
    std::map<int, double> res;
    const double* Mz;/// Z-boson mass
    const double* alphas;/// alphaS(Mz^2)
  };

extern "C" {
    // Fortran routines
#define pbtmdsubr pbtmdsubr_
    double pbtmdsubr(int &i, double &x, double  &q2, double*  mc, double*  mb, char* name, int &len);
#define write_pbtmd write_pbtmd_
    double write_pbtmd(const char* name, int &len);
    
    }
};
