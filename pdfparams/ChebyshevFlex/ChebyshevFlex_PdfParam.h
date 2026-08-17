#pragma once

#include "BasePdfParam.h"

/**
  @class ChebyshevFlex_PdfParam

  @brief A highly optimized, flexible Chebyshev pdf parameterisation

  Формула:
  xf(x) = p1 * x^p2 * (1-x)^p3 * [ 1 + Sum_{i=1}^{n-4} (p_{i+3} * T_i(y(x))) ]
  
  Мапування:
  y(x) = 1 - 2 * x^k
  
  де k = p0 (нульовий параметр).
*/

namespace xfitter {
class ChebyshevFlex_PdfParam : public BasePdfParam {
  public:
    ChebyshevFlex_PdfParam(const std::string& name) : BasePdfParam(name) {}
    virtual double operator()(double x) const override final;
    virtual void atStart() override final;
};
}