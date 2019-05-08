
#pragma once

#include "BasePdfParam.h"

/**
  @class ABMPseaPdfParam

  @brief A class for ABMPsea pdf parameterisation according to Eqs. 19-22 from Phys.Rev. D96 (2017) no.1, 014011
  xv(x) = A * (1 - x)^b * x^[a * (1 + gam_m1 * ln(x)) * (1 + gam_1 * x + gam_2 * x^2 + gam_3 * x^3)]

  @version 0.1
  @date 2019-02-25
  */

namespace xfitter{
  class ABMPseaPdfParam:public BasePdfParam{
  public:
    ABMPseaPdfParam(const std::string&inName):BasePdfParam(inName){}
    //Evaluate xf(x) at given x with current parameters
    virtual double operator()(double x)const override final;
    // (Optional) compute moments:
    // virtual double  moment(int nMoment=-1)const override final;
    // (Optional) set moments:
    // virtual void setMoment(int nMoment,double value)override final;
    // (Optional)
    //Initialize from a yaml node. Uses node[getName] as the basis
    // virtual void initFromYaml(YAML::Node value)override final;
  };
}
