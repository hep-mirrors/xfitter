
#pragma once

#include "BasePdfParam.h"

/**
  @class Fantomas_PdfParam 

  @brief A class for Fantomas_ pdf parameterisation

  HERAPDF - style parameterisation:
  xf(x)=A*x^B*(1-x)^C*(1+P(x))
  where P(x) is a polynomial with other parameters as coefficients

  Number of parameters may vary, but at least 3, corresponding to A,B,C

  In terms of indexed parameters:
  xf(x)=par_0*x^par_1*(1-x)^par_2*(1+sum_{i=3}^{N}{par_i*x^(i-2)})

  @version 0.2
  @date 2018-08-14
  */

namespace xfitter{
class Fantomas_PdfParam:public BasePdfParam{
  public:
    Fantomas_PdfParam(const std::string&inName):BasePdfParam(inName){}
    virtual double operator()(double x)const override final;
    virtual double moment(int nMoment=-1)const override final;
    virtual void setMoment(int nMoment,double value)override final;
    virtual void atStart()override final;
    virtual void atIteration()override final;
 private:
    void updateParameters();
};
}
