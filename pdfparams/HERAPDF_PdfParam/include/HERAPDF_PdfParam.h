
#pragma once

#include "BasePdfParam.h"

/**
  @class HERAPDF_PdfParam 

  @brief A class for HERAPDF_ pdf parameterisation

  Based on the ReactionTheory class. Reads options produces 3d cross section.

  @version 0.1
  @date 2018-07-11
  */

class HERAPDF_PdfParam : public BasePdfParam
{
  public:
     /// Default constructor. Name is the PDF name
     HERAPDF_PdfParam (const std::string& inName) : BasePdfParam(inName) {}
     /// Compute xf(x,pars). Pure virtual method
     virtual double compute( double const x, double const* pars) const override final;
    
     /// (optionally) compute moments:
     // virtual double moment( double const* pars, int const iMoment = 1) override final;
};
