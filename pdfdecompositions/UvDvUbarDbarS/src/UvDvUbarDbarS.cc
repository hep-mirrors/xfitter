 
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
#include <iomanip>
#include <cmath>

/// TEMPORARY XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
extern "C" {
  /// Interface to minuit parameters
  void addexternalparam_(const char name[],  const double &val, 
                         const double  &step,
                         const double &min, const double &max, 
                         const double &prior, const double &priorUnc,
                         const int &add, 
                         map<std::string,double*> *map,
                         int len);
}


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


      /// HARDWIRE old-way for now:  XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      for (int i = 0; i< pParam->getNPar(); i++) {
	// get a good name
	const std::string pnam = getName() + "_"+ pdfName + "_p" +  std::to_string(i) ;
	std::cout << pnam ;
	std::cout << " Pars: " << parValues[i] <<std::endl;

	double val = parValues[i];
	double step = std::fabs(val)/100.;       /// if 0, parameter is fixed !!! 

	
	double minv  = 0;
	double maxv  = 0;
	double priorVal = 0;
	double priorUnc = 0;
	int add = true;

	
	/// Here it goes to minuit:
	addexternalparam_(pnam.c_str(), val, step, minv, maxv, priorVal, priorUnc, add, &XFITTER_PARS::gParameters,  pnam.size() );
      }
      
      addParameterisation(pdfName,pParam);
    }
  }

  void UvDvUbarDbarS::initAtIteration() {
    // counting sum-rules for uv and dv
    const BasePdfParam* pdfParam = getPdfParam("xuv");
    std::unique_ptr<double[]> pars = getParValues(pdfParam);  
    pars[0] = 1.;
    double sum0 = pdfParam->moment(pars.get(),-1);
    _uSum = 2.0 / sum0;

    pdfParam = getPdfParam("xdv");
    pars = getParValues(pdfParam);  
    pars[0] = 1.;
    sum0 = pdfParam->moment(pars.get(),-1);
    _dSum = 1.0 / sum0;

    // momentum sum-rule

    // quark part
    double xsumq = 0;
    const std::vector<std::string> nameS{"xubar","xdbar","xs","xuv","xdv"};
    for ( auto const& name : nameS ) {
      const BasePdfParam* pdfParam = getPdfParam(name);
      std::unique_ptr<double[]> pars = getParValues(pdfParam);

      if ( name == "xuv") {
	pars[0] = _uSum;
      }
      if (name == "xdv" ) {
	pars[0] = _dSum;
      }
      xsumq += pdfParam->moment(pars.get(),0);
    }
    
    // gluon part
    pdfParam = getPdfParam("xg");
    pars = getParValues(pdfParam);
    pars[0] = 1.;
    double xsumg = pdfParam->moment(pars.get(),0);
    
    _gSum = (1. - xsumq)/xsumg;

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

      if (name == "xg") {
	pars[0] = _gSum;
      }
      if (name == "xuv") {
	pars[0] = _uSum;
      }
      if (name == "xdv") {
	pars[0] = _dSum;
      }
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
    for (int i=0; i<param->getNPar(); i++) {
      const std::string pnam = getName() + "_" + param->getName() + "_p" +  std::to_string(i) ;
      pars[i] = *XFITTER_PARS::gParameters[pnam];
    }
    return pars;
  }
  
  
  double UvDvUbarDbarS::valence(double x, std::string const& name, double sum) const {
    const BasePdfParam* pdfParam = getPdfParam(name);
    // get parameters
    std::unique_ptr<double[]> pars = getParValues(pdfParam);  

    pars[0] = sum;
    
    return pdfParam->compute(x,pars.get());
  }

  double UvDvUbarDbarS::sea(double x, std::string const& name) const {
    const BasePdfParam* pdfParam = getPdfParam(name);
    // get parameters
    std::unique_ptr<double[]> pars = getParValues(pdfParam);  
    
    return pdfParam->compute(x,pars.get());
  }

  
  double UvDvUbarDbarS::uv(double x) const {
    return valence(x,"xuv",_uSum); 
  }

  double UvDvUbarDbarS::dv(double x) const {
    return valence(x,"xdv",_dSum); 
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

  double UvDvUbarDbarS::g(double x) const {
    const BasePdfParam* pdfParam = getPdfParam("xg");
    // get parameters
    std::unique_ptr<double[]> pars = getParValues(pdfParam);  
    pars[0] = _gSum;
    return pdfParam->compute(x,pars.get());
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
	{-2,dbar(x)},
	{-1,ubar(x)},
	{ 1,ubar(x)+uv(x)},
	{ 2,dbar(x)+dv(x)},
	{ 3,s(x)},
	{ 4,0.},
	{ 5,0.},
	{ 6,0.},
	{21,g(x)}
      };
      return res;
    };
    return _f0;
  }
}


