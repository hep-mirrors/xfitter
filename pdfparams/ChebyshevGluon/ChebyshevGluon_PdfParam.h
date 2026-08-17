#pragma once

#include "BasePdfParam.h"

namespace xfitter
{
  class ChebyshevGluon : public BasePdfParam
  {
  public:
    ChebyshevGluon(const char* name);
    
    // Прибрали 'override final', бо цього методу немає в BasePdfParam
    virtual const char* getClassName() const;
    
    virtual double operator()(double x) const override final;
  };
}