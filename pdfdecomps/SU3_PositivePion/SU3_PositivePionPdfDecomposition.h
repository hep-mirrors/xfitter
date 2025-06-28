#pragma once
#include "BasePdfDecomposition.h"
#include"BasePdfParam.h"
/**
  @class SU3_PositivePionPdfDecomposition 

  @brief A class for pdf decomposition for the pion with SU3-symmetric sea

  Used for pi-
  Assumes that at starting scale:
    ubar=d
    dbar=u=s=sbar
  Parametrised distributions are:
    v:=dval-uval=2(d-u)
    S:=2u+2dbar+s+sbar=6u
    g
  Therefore, transformations to flavor basis:
    d=ubar=v/2+S/6
    u=dbar=s=sbar=S/6
    g=g
    others=0
  And sum rules for pi- are:
    \int_0^1 v dx=2
    \int_0^1 x*(v+S+g) dx=1
  @version 0.2
  @date 2018-08-14
  */
namespace xfitter{
class SU3_PositivePionPdfDecomposition:public BasePdfDecomposition{
  public:
    SU3_PositivePionPdfDecomposition(const char*name);
    virtual const char*getClassName()const override final;
    virtual void atStart()override final;
    virtual void atIteration()override final;
    virtual std::map<int,double>xfxMap(double x)const override final;
  private:
    BasePdfParam*par_v{nullptr},*par_S{nullptr},*par_g{nullptr};
  };
}
