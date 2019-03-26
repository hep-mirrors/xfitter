
#pragma once

#include "BasePdfParam.h"

/**
  @class PolySqrtPdfParam 

  @brief A class for PolySqrt pdf parameterisation

  xf(x)=A*x^B*(1-x)^C*(1+P(sqrt(x)))
  where P(sqrt(x)) is a polynomial of sqrt(x) with other parameters as coefficients:
  P(sqrt(x))=D*sqrt(x)+E*x+F*x*sqrt(x)+...

  Number of parameters may vary, but at least 3, corresponding to A,B,C

  In terms of indexed parameters:
  xf(x)=par_0*x^par_1*(1-x)^par_2*(1+sum_{i=3}^{N}{par_i*x^(i/2-1)})

  @version 0.1
  @date 2018-08-16
  */

namespace xfitter{
class PolySqrtPdfParam:public BasePdfParam{
  public:
    PolySqrtPdfParam(const std::string&inName):BasePdfParam(inName){}
    virtual double operator()(double x)const override final;
    virtual double moment(int nMoment=-1)const override final;
    virtual void atStart()override final;
};
}
