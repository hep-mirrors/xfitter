#pragma once
#include"BasePdfParam.h"
#include"tinyexpr.h"

/**
  @class ExpressionPdfParam

  @brief A class for Expression pdf parameterisation

  Write parameterisation as an expression directly in input YAML file.
  Based on TinyExpr parser.
  Expression may contain math functions supported by TinyExpr (sin, log, exp, ...)
  , minimization parameters and 'x'
  Note that '^' is exponentiation.

  For sumrules, numerical integration as implemented in BasePdfParam is used

  If given expression has form "PAR*rest_of_expression",
  parameter PAR will be scaled to enforce sumrules (to set moment)
  Otherwise moment of expression cannot be set, and trying to do so will raise an error

  Example:
  Parameterisations:
    example:
      class: Expression
      expression: "Av*x^Bv*(1-x)^Cv"
*/

namespace xfitter{
class ExpressionPdfParam:public BasePdfParam{
  public:
    ExpressionPdfParam(const std::string&name):BasePdfParam(name){}
    virtual double operator()(double x)const override final;
    virtual void setMoment(int nMoment,double value)override final;
    virtual void atStart()override final;
  private:
    te_expr*expr=nullptr;
    double*normPar=nullptr;//normalization parameter, will be modified to set moment
    //normPar==nullptr <=> cannot set moment
    mutable double x;
};
}
