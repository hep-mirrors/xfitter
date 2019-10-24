#pragma once
#include"BaseEvolution.h"

namespace xfitter{
/**
  @class FlipUD_Evol

  @brief A class for FlipUD evolution

  Takes as input another evolution and switches up-quark with down-quark, ubar with dbar
  This turns a proton into neutron and vice versa, assuming isospin symmetry

  In YAML configuration, provide another evolution's name as "input"
*/
class FlipUD_Evol:public BaseEvolution{
  public:
    FlipUD_Evol(const char*name):BaseEvolution(name){}
    virtual const char*getClassName()const override final{return"FlipUD";};
    virtual void atStart()override final;
    virtual void atConfigurationChange()override final;
    virtual std::map<int,double>xfxQmap(double x,double Q)override final;
    virtual double xfxQ(int i,double x,double Q)override final;
    virtual void xfxQarray(double x,double Q,double*pdfs)override final;
    virtual double getAlphaS(double Q)override final;
    virtual std::vector<double> getXgrid() override final;
    virtual std::vector<double> getQgrid() override final;
  private:
    BaseEvolution*input;
};
}
