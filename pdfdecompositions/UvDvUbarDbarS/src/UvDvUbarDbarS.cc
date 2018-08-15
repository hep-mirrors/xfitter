 
/*
   @file UvDvUbarDbarS_PdfDecomposition.cc
   @date 2018-08-14
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-11
*/

#include "UvDvUbarDbarS.h"
#include "HERAPDF_PdfParam.h"
#include "xfitter_pars.h"
#include <iostream>
#include <iomanip>
#include <cmath>

namespace xfitter
{
  /// the class factories
  extern "C" UvDvUbarDbarS* create() {
    return new UvDvUbarDbarS();
  }

  
  //_________________________________________________________________________________
  UvDvUbarDbarS::UvDvUbarDbarS(): BasePdfDecomposition{"UvDvUbarDbarS"} { }

  //_________________________________________________________________________________
  void UvDvUbarDbarS::initAtStart(const std::string & pars){
    // HARDWIRE A BIT FOR NOW XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    // Currently we create and initialise parameterisations
    // This is not great, it is not Decomposition's job, they should be created globally
    const YAML::Node pdfs    = XFITTER_PARS::gParametersY.at("UvDvubardbars")["HERAPDF_pdfparam"];

    for (const auto& pdf : pdfs) {
      const string pdfName = pdf.first.as<string>();
      BasePdfParam* pParam = new HERAPDF_PdfParam(pdfName);
      pParam->initFromYaml(pdf.second);
			addParameterisation(pdfName,pParam);
    }
  }

  void UvDvUbarDbarS::initAtIteration() {
    //Enforce sum rules
    // counting sum-rules for uv and dv
    getPdfParam("xuv")->setMoment(-1,2.0);
    getPdfParam("xdv")->setMoment(-1,1.0);
    // momentum sum-rule
    // quark part
    double xsumq=0;
    xsumq+=  getPdfParam("xuv"  )->moment(0);
    xsumq+=  getPdfParam("xdv"  )->moment(0);
    xsumq+=2*getPdfParam("xubar")->moment(0);
    xsumq+=2*getPdfParam("xdbar")->moment(0);
    xsumq+=2*getPdfParam("xs"   )->moment(0);
    // gluon part
    getPdfParam("xg")->setMoment(0,xsumq);

    printParams();
  }

  void UvDvUbarDbarS::printParams() {
    std::cout << "\n" << std::left<< std::setw(8) << " Name " << std::right;
    for ( int i =0 ; i<6; i++) {
      std::cout << std::setw(12) << " par"+std::to_string(i) ;
    }
    std::cout << "\n";
    for ( auto p : _pdfParams) {
      auto pdfParam = p.second;
      auto name = p.first;
      auto pars = getParValues(pdfParam);
      int  npar = pdfParam->getNPar();

      //      std::cout << "name :" << name << " npar: " << npar << "\n";
      std::cout << std::left<< std::setw(8) << name << std::right;
      for ( int i =0 ; i<npar; i++) {
        std::cout << std::setw(12) << pars[i];
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  
  std::unique_ptr<double[]> UvDvUbarDbarS::getParValues(BasePdfParam const* param) const {
    std::unique_ptr<double[]> pars( new double[param->getNPar()] );
    for (unsigned int i=0; i<param->getNPar(); i++) {
      const std::string pnam=param->getName()+ "_p" +  std::to_string(i) ;
      pars[i] = *XFITTER_PARS::gParameters.at(pnam);
    }
    return pars;
  }
  
  double UvDvUbarDbarS::valence(double x, const std::string&name) const {
    return(*getPdfParam(name))(x);
  }

  double UvDvUbarDbarS::sea(double x,const std::string&name)const{
    return(*getPdfParam(name))(x);
  }

  double UvDvUbarDbarS::uv(double x) const {
    return valence(x,"xuv"); 
  }

  double UvDvUbarDbarS::dv(double x) const {
    return valence(x,"xdv"); 
  }

  double UvDvUbarDbarS::ubar(double x) const {
    return sea(x,"xubar");
  }

  double UvDvUbarDbarS::dbar(double x) const {
    return sea(x,"xdbar"); // XXXXXXXXXXXXX HARDWIRE
  }
  
  double UvDvUbarDbarS::s(double x) const {
    return sea(x,"xs");
  }

  double UvDvUbarDbarS::g(double x)const{
    return (*getPdfParam("xg"))(x);
  }

  
  //_________________________________________________________________________________
  std::function<std::map<int,double>(const double& x)> UvDvUbarDbarS::f0() const
  {
    // lambda function
    const auto _f0 = [=] (double const& x)->std::map<int, double> {
      std::map<int, double> res  = {
        {-6,0},
        {-5,0},
        {-4,0},
        {-3,s(x)},
        {-2,ubar(x)},
        {-1,dbar(x)},
        { 1,dbar(x)+dv(x)},
        { 2,ubar(x)+uv(x)},
        { 3,s(x)},
        { 4,0},
        { 5,0},
        { 6,0},
        {21,g(x)}
      };
      return res;
    };
    return _f0;
  }
}


