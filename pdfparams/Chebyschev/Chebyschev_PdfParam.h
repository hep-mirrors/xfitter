
#pragma once

#include "BasePdfParam.h"

/**
  @class Chebyschev_PdfParam 

  @brief A class for Chebyschev_ pdf parameterisation

  Chebyschev - style parameterisation:
  xf(x)=A*x^B*(1-x)^C*(1+P(y(x)))
  where P(y(x)) = 1 + Sum ( a_i*T_i(y(x)) )
  with Tch_i the Chebyschev polynomial of order i,
  y(x) = 1 - 2*x^k with k=0.5
  and the a_i parameters as coefficients

  Number of parameters may vary, but at least 3, corresponding to A,B,C

  In terms of indexed parameters:
  xf(x)=par_0*x^par_1*(1-x)^par_2*(1+sum_{i=3}^{N}{par_i*x^(i-2)})

  @version 0.1
  @date 2023-01-14
  */

namespace xfitter{
  class Chebyschev_PdfParam:public BasePdfParam{
  public:
    Chebyschev_PdfParam(const std::string&inName):BasePdfParam(inName){}
    virtual double operator()(double x) const override final;
    virtual void atStart() override final;
  private:
    double Tn(unsigned int n, double x) const;
  };
}
