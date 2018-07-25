 
/*
   @file UvDvUbarDbarS_PdfDecomposition.cc
   @date 2018-07-11
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-11
*/

#include "UvDvUbarDbarS.h"
#include "HERAPDF_PdfParam.h"
#include "xfitter_pars.h"
#include <iostream>

namespace xfitter
{
  /// the class factories
  extern "C" UvDvUbarDbarS* create() {
    return new UvDvUbarDbarS();
  }

  
  //_________________________________________________________________________________
  UvDvUbarDbarS::UvDvUbarDbarS(): BasePdfDecomposition{"UvDvUbarDbarS"} { }

  //_________________________________________________________________________________
  void UvDvUbarDbarS::initAtStart(const std::string & pars)
  {
    // HARDWIRE A BIT FOR NOW XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    const YAML::Node pdfs    = XFITTER_PARS::gParametersY.at("UvDvubardbars")["HERAPDF_pdfparam"];

    for (const auto& pdf : pdfs) {
      const string pdfName = pdf.first.as<string>();
      BasePdfParam* pParam = new HERAPDF_PdfParam(pdfName);
      // read its parameters:
      double *parValues = pParam->initFromYaml(pdf.second);
      std::cout << " Pars: " << parValues[0] <<std::endl;
      addParameterisation(pParam);
    }
  }

  double UvDvUbarDbarS::uv(double x) const {
    return 0;
  }

  double UvDvUbarDbarS::dv(double x) const {
    return 0;
  }

  double UvDvUbarDbarS::dbar(double x) const {
    return 0;
  }

  double UvDvUbarDbarS::ubar(double x) const {
    return 0;
  }
  
  double UvDvUbarDbarS::s(double x) const {
    return 0;
  }

  
  //_________________________________________________________________________________
  std::function<std::map<int,double>(const double& x)> UvDvUbarDbarS::f0() const
  {
    // lambda function
    const auto _f0 = [=] (double const& x)->std::map<int, double> {
      double ubar_ = ubar(x);
      double dbar_ = dbar(x);
      double u_ = ubar_+uv(x);
      double d_ = dbar_+dv(x);
      double s_ = s(x);
      std::map<int, double> res  = {
	{-6,0},	
	{-5,0},
	{-4,0},
	{-3,s_},
	{-2,dbar_},
	{-1,ubar_},
	{ 1,u_},
	{ 2,d_},
	{ 3,s_},
	{ 4,0.},
	{ 5,0.},
	{ 6,0.}
      };
      return res;
    };
    return _f0;
  }
}


