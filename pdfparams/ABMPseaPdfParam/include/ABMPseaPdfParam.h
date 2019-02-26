
#pragma once

#include "BasePdfParam.h"

/**
  @class ABMPseaPdfParam 

  @brief A class for ABMPsea pdf parameterisation

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