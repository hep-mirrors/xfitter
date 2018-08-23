 
/*
   @file EvolutionLHAPDF.cc
   @date  2018-08-20
   @author  AddEvolution.py
   Created by  AddEvolution.py on 2018-08-20
*/

#include "EvolutionLHAPDF.h"
#include "xfitter_pars.h"
#include "CheckForPDF.h"

namespace xfitter
{

// the class factories
extern "C" EvolutionLHAPDF* create() {
  return new EvolutionLHAPDF();
}

    
/// Global initialization
  void EvolutionLHAPDF::initAtStart() {

    return ;
 };
    
  /// Init at each iteration
  void EvolutionLHAPDF::initAtIteration() {

    // Initialize LHAPDF
    auto pars = XFITTER_PARS::gParametersY["LHAPDF"];

    _set_name = pars["set"].as<std::string>();
    _member   = pars["member"].as<int>();

    // check if exists first
    CheckForPDF(_set_name.c_str());
    
    _pdf      = LHAPDF::mkPDF(_set_name,_member);
    return ;
  };

  /// Return PDFs as a map <int,double> where int is PDF ID (-6, ... 6, 21)   
  std::function<std::map<int,double>(double const& x, double const& Q)> EvolutionLHAPDF::xfxQMap() {
    auto f0 = [=] (double const& x, double const& Q)->std::map<int, double>{ return _pdf->xfxQ(x, Q); };
    return  f0; 
  };

  /// Returns PDFs as a function of i, x, Q
  std::function<double(int const& i, double const& x, double const& Q)> EvolutionLHAPDF::xfxQDouble() {
    auto f0 = [=] (int const& i, double const& x, double const& Q)->double { return _pdf->xfxQ(i,x, Q); };
    return  f0; 
  };
    
  /// Returns PDFs as double pdfs* --> double[13] from -6 to 6.  
  std::function<void(double const& x, double const& Q, double* pdfs)> EvolutionLHAPDF::xfxQArray() {
    auto f0 = [=] (double const& x, double const& Q, double* pdfs)->void
      {
	std::vector<double> vpdfs{13};	
	_pdf->xfxQ(x, Q, vpdfs);
	for (int i=0; i<13; i++) 
	  pdfs[i] = vpdfs[i];
	return ;
      };
    return  f0; 
  };

  /// Returns alphaS
  std::function<double(double const& Q)> EvolutionLHAPDF::AlphaQCD() {
    auto f0 = [=] ( double const& Q) -> double {
      return _pdf->alphasQ(Q);
    };
    return f0;
  };
}
 
