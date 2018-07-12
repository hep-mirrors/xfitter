
#pragma once

#include "BasePdfDecomposition.h"
#include "BasePdfParam.h"

/**
  @class UvDvUbarDbarS_PdfDecomposition 

  @brief A class for UvDvUbarDbarS_ pdf decomposition

  @version 0.1
  @date 2018-07-11
  */

class UvDvUbarDbarS_PdfDecomposition : public BasePdfDecomposition
{
 public:
     /// Default constructor. Name is the PDF name
    UvDvUbarDbarS_PdfDecomposition (const std::string& inName);

    /// Optional initialization at the first call
    virtual void initAtStart(const std::string & pars) override final;

     /// Compute PDF in a physical base in LHAPDF format for given x and Q
    virtual std::map<int,double> compute ( const double& x, const double& Q) const override final; 

 private:
    std::map<std::string,BasePdfParam*> pdf_pars;
};
