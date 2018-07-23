 
/*
   @file UvDvUbarDbarS_PdfDecomposition.cc
   @date 2018-07-11
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-11
*/

#include "UvDvUbarDbarS.h"
#include "HERAPDF_PdfParam.h"
#include "xfitter_pars.h"

namespace xfitter
{
  //_________________________________________________________________________________
  UvDvUbarDbarS::UvDvUbarDbarS(): BasePdfDecomposition{"UvDvUbarDbarS"} { }

  //_________________________________________________________________________________
  void UvDvUbarDbarS::initAtStart(const std::string & pars)
  {
    // HARDWIRE A BIT FOR NOW XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    const YAML::Node pdfs    = XFITTER_PARS::gParametersY.at("UvDvUbarDbarS")["HERAPDF_pdfparam"];
    for (const auto& pdf : pdfs) {
      const string pdfName = pdf.first.as<string>();
      BasePdfParam* pParam = new HERAPDF_PdfParam(pdfName);
      // read its parameters:
      double *parValues;
      pParam->initFromYaml(pdf,parValues);
      addParameterisation(pParam);
    }
  }

  //_________________________________________________________________________________
  std::function<std::map<int,double>(const double& x)> UvDvUbarDbarS::f0() const
  {
    return nullptr;
  }
}


