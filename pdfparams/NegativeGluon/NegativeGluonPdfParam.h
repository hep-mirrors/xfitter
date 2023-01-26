
#pragma once

#include "BasePdfParam.h"

/**
  @class NegativeGluonPdfParam 

  @brief A class for NegativeGluonPdfParam pdf parameterisation. Uses HERAPDF-style 5 parameters for positive and 3 for negative gluon

  HERAPDF - style parameterisation:
  xf(x)= A x^B (1-x)^C (1+D*x+E*x^2) - A' x^B' (1-x)^C'

  Need exactly 8 parameters, corresponding to A,B,C,D,E,A',B',C'

  As defined for HERAPDF, when performing the sum rule calculation, A is changed, but A' is left unchanged

  @version 0.2
  @date 2022-12-09
  */

namespace xfitter {
  class NegativeGluonPdfParam:public BasePdfParam{
  public:
    NegativeGluonPdfParam(const std::string&inName):BasePdfParam(inName){}
    //Evaluate xf(x) at given x with current parameters
    virtual double operator()(double x)const override final;
    virtual double moment(int nMoment=-1)const override final;
    virtual void atStart()override final;
    virtual void setMoment(int nMoment,double value) override;

  private:
    double moment_pos(int nMoment=-1) const;
    double moment_neg(int nMoment=-1) const;
  };
};
