
#pragma once

#include "BasePdfParam.h"

/**
  @class NegativeGluonPdfParam 

  @brief A class for NegativeGluonPdfParam pdf parameterisation. Uses HERAPDF-style 5 parameters for positive and 5 for negative gluon

  @version 0.1
  @date 2018-09-26
  */

namespace xfitter {
  class NegativeGluonPdfParam:public BasePdfParam{
  public:
    NegativeGluonPdfParam(const std::string&inName):BasePdfParam(inName){}
    //Evaluate xf(x) at given x with current parameters
    virtual double operator()(double x)const override final;
    virtual double moment(int nMoment=-1)const override final;
    virtual void atStart()override final;
  };
};
