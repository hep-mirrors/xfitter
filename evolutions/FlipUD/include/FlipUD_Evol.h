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
    virtual std::function<std::map<int,double>(double const&x,double const&Q)>xfxQMap()override final;
    virtual std::function<void(double const&x,double const&Q,double*pdfs)>xfxQArray()override final;
    virtual std::function<double(int const&i,double const&x,double const&Q)>xfxQDouble()override final;
    virtual std::function<double(double const&Q)>AlphaQCD()override final;
  private:
    BaseEvolution*input;
};
}
