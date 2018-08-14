
#pragma once

#include "BasePdfParam.h"

/**
  @class HERAPDF_PdfParam 

  @brief A class for HERAPDF_ pdf parameterisation

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.2
  @date 2018-08-14
  */

class HERAPDF_PdfParam:public BasePdfParam{
  public:
    HERAPDF_PdfParam(const std::string&inName):BasePdfParam(inName){}
    virtual double operator()(double x)const override final;
};
