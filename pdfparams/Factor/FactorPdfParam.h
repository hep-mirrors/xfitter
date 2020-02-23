#pragma once
#include"BasePdfParam.h"

/**
  @class FactorPdfParam

  @brief A class for Factor pdf parameterisation

  Takes another parameterisation as input and multiplies it by a factor, which is provided as a fittable parameter

  YAML settings:
  input: name of another parameterisation
  factor: name of parameter to use as factor

*/

namespace xfitter{
class FactorPdfParam:public BasePdfParam{
  public:
    FactorPdfParam(const std::string&name):BasePdfParam(name){}
    virtual double operator()(double x)const override final;
    virtual double  moment(int nMoment=-1)const override final;
    virtual void setMoment(int nMoment,double value)override final;
    virtual void atStart()override final;
  private:
    BasePdfParam* input;
    const double* factor;
};
}
