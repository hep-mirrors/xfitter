#include "EvolutionHoppet.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include <iostream>

#include "hoppet_v1.h"

namespace xfitter
{
  // the class factories
  extern "C" EvolutionHoppet*create(const char*name){
    return new EvolutionHoppet(name);
  }
  const char*EvolutionHoppet::getClassName()const{return "Hoppet";}


  //_________________________________________________________________________________
  void EvolutionHoppet::atStart()
  {
    std::cout << " HERE WE ARE in HOPPET " << std::endl;
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    _inPDFs=XFITTER_PARS::getInputDecomposition(yamlNode);
    // Retrieve parameters needed to initialize APFEL++.
    const double* MCharm   = XFITTER_PARS::getParamD("mch");
    const double* MBottom  = XFITTER_PARS::getParamD("mbt");
    const double* MTop     = XFITTER_PARS::getParamD("mtp");
    const YAML::Node xGrid = yamlNode["xGrid"];

    // Add hoppetStart
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
    double dy = 0.1;
    hoppetStart(dy, PtOrder);
    
    atConfigurationChange();
  }

  //_________________________________________________________________________________
  void EvolutionHoppet::atIteration()
  {
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    // Retrieve the relevant parameters needed to compute the evolutions
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;
    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    const double* Q_ref      = XFITTER_PARS::getParamD("Mz");
    const double* Alphas_ref = XFITTER_PARS::getParamD("alphas");
    const YAML::Node QGrid   = yamlNode["QGrid"];

    // add hoppetEvolve
  }
  std::map<int,double>EvolutionHoppet::xfxQmap(double x,double Q){
  }

  double EvolutionHoppet::xfxQ(int i,double x,double Q){
  }

  void EvolutionHoppet::xfxQarray(double x,double Q,double*pdfs){
  }

  double EvolutionHoppet::getAlphaS(double Q){
  }

  vector<double> EvolutionHoppet::getXgrid() {
  }

  vector<double> EvolutionHoppet::getQgrid() {
  }
}
