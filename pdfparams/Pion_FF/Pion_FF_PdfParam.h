
#pragma once

#include "BasePdfParam.h"

/**
  @class Pion_FF_PdfParam 

  @brief A class for Pion-FF pdf parameterisation
  @version 0.2
  @date 2018-08-14
  */

namespace xfitter{
class Pion_FF_PdfParam:public BasePdfParam{
 public:
  Pion_FF_PdfParam(const std::string&inName):BasePdfParam(inName){}
  virtual double operator()(double x)const override final;
  virtual void atStart()override final;
  virtual void atIteration()override final;
};
}
