#include "EvolutionHOPPET.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include <iostream>
#include <cmath>
#include "hoppet_v1.h"

using namespace std;
using namespace hoppetv1;

// Global var to hold current pdfDecomposition
xfitter::BasePdfDecomposition *gPdfDecomp = nullptr;

namespace xfitter
{

  // the class factories
  extern "C" EvolutionHOPPET*create(const char*name){
    return new EvolutionHOPPET(name);
  }
  const char*EvolutionHOPPET::getClassName()const{return "HOPPET";}

  void EvolutionHOPPET::atStart()
  {
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    _inPDFs=XFITTER_PARS::getInputDecomposition(yamlNode);
    //const YAML::Node xGrid = yamlNode["xGrid"];
    
    _order = OrderMap(XFITTER_PARS::getParamS("Order"));
    // temporary: allow different orders in evolution and DIS SFs
    if (XFITTER_PARS::gParametersS.find("Order_HOPPET_Evolution") != XFITTER_PARS::gParametersS.end()) {
      _order = OrderMap(XFITTER_PARS::getParamS("Order_HOPPET_Evolution"));
    }

    double dy = yamlNode["dy"].as<double>();
    //hoppetStart(dy, _order);
    double dy_over_dlnlnQ = yamlNode["dy_over_dlnlnQ"].as<double>();
    double dlnlnQ = dy/dy_over_dlnlnQ;
    _ymax = yamlNode["ymax"].as<double>();
    _Qmin = yamlNode["Qmin"].as<double>();
    _Qmax = yamlNode["Qmax"].as<double>();
    int order_interpol = yamlNode["order_interpol"].as<int>();
    hoppetStartExtended(_ymax, dy, _Qmin, _Qmax, dlnlnQ, _order, order_interpol, factscheme_MSbar);
    _isFFNS = 0; // VFNS by default
    if (yamlNode["isFFNS"]) {
      _isFFNS = yamlNode["isFFNS"].as<int>();
    }
    else {
      if(XFITTER_PARS::gParametersI.find("isFFNS") != XFITTER_PARS::gParametersI.end()) {
        _isFFNS = XFITTER_PARS::gParametersI.at("isFFNS");
      }
    }
    _msbar = 0; // pole mass scheme by default
    if (yamlNode["MSbar"]) {
      _msbar = yamlNode["MSbar"].as<int>();
    }
    else {
      if(XFITTER_PARS::gParametersI.find("MSbar") != XFITTER_PARS::gParametersI.end()) {
        _msbar = XFITTER_PARS::gParametersI.at("MSbar");
      }
    }
    _nflavour = -1;
    if (yamlNode["NFlavour"]) {
      _nflavour = yamlNode["NFlavour"].as<int>();
    }
    else {
      if(XFITTER_PARS::gParametersI.find("NFlavour") != XFITTER_PARS::gParametersI.end()) {
        _nflavour = XFITTER_PARS::gParametersI.at("NFlavour");
      }
    }
    if(_isFFNS == 1) {
      hoppetSetFFN(_nflavour);
    }
    else if(_isFFNS == 0) {
      _mch = XFITTER_PARS::getParamD("mch");
      _mbt = XFITTER_PARS::getParamD("mbt");
      _mtp = XFITTER_PARS::getParamD("mtp");
      //hoppetSetVFN(*_mch, *_mbt, *_mtp);
      hoppetSetPoleMassVFN(*_mch, *_mbt, *_mtp);
      if (_msbar) {
        hoppetSetMSbarMassVFN(*_mch, *_mbt, *_mtp);
      }
    }
    else
    {
      hf_errlog(20240903, "F: Unsupported isFFNS = " + std::to_string(_isFFNS));
    }
    atConfigurationChange();
    printf("SZ HOPPET _isFFNS %d MSbar %d _nflavour %d\n", _isFFNS, _msbar, _nflavour);
  }

void  heralhc_init(const double & x,
                    const double & Q,
                    double * pdf)                  
 {
  const std::map<int, int> ip =
      {{-5, -5},{-4, -4},{-3, -3}, {-2, -2}, {-1, -1}, {0, 21}, {1, 1}, {2, 2}, {3, 3},{4, 4},{5, 5}};
  for (auto &i : ip)
  {
    pdf[i.first + 6] = gPdfDecomp->xfxMap(x)[i.second];
    //double t = gPdfDecomp->xfxMap(x)[i.second]; 
    //std::cout << " from decomposition: " << i.first << " " << i.second << "  " << t << std::endl;
  }
 } 
  
    void EvolutionHOPPET::atIteration()
  {
    // use  https://github.com/hoppet-code/hoppet/blob/master/example_f77/cpp_tabulation_example.cc
    //std::cout << " HERE WE ARE in HOPPET " << std::endl;
    //exit(0);
    
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    // Retrieve the relevant parameters needed to compute the evolutions
    double* Q0 = XFITTER_PARS::getParamD((yamlNode["Q0"]) ? yamlNode["Q0"].as<string>() : "Q0");
    try {
      _alphas_q0 = XFITTER_PARS::getParamD((yamlNode["alphas_Q0"]) ? yamlNode["alphas_Q0"].as<string>() : "alphas_Q0");
    }
    catch(std::out_of_range&ex) {
      _alphas_q0 = nullptr;
    }
    if (!_alphas_q0) {
      _alphas_q0 = XFITTER_PARS::getParamD("Mz");
    }
    _alphas = XFITTER_PARS::getParamD((yamlNode["alphas"]) ? yamlNode["alphas"].as<string>() : "alphas");
    //const YAML::Node QGrid   = yamlNode["QGrid"];
    gPdfDecomp = _inPDFs;

    if(_isFFNS == 0) {
      //hoppetSetVFN(*_mch, *_mbt, *_mtp);
      hoppetSetPoleMassVFN(*_mch, *_mbt, *_mtp);
      if (_msbar) {
        hoppetSetMSbarMassVFN(*_mch, *_mbt, *_mtp);
      }
    }

    hoppetEvolve( *_alphas, *_alphas_q0, _order, 1.0, heralhc_init, *Q0);
    //double f[13];
    //hoppetEval(0.001,10.,f);
    //printf("f[6] = %f  as = %f\n", f[6], hoppetAlphaS(10.));
    //std::cout << " HERE WE ARE OUT OF HOPPET " << std::endl;
  }
  
  
  std::map<int,double>EvolutionHOPPET::xfxQmap(double x,double Q){
    //std::cout << " HERE WE ARE in HOPPET A " << std::endl;
    double pdfs[14];
    xfxQarray(x, Q, pdfs);
    std::map<int, double> res;
    
    const int npdfMax = 6;

    for (int ipdf = -6; ipdf <= npdfMax; ipdf++)
      {
	int ii = (ipdf == 0) ? 21 : ipdf;
	// photon PDF:
	if (ipdf == 7)
	  ii = 22;
	res[ii] = pdfs[ipdf+6];
      }
    return res;
  }

  double EvolutionHOPPET::xfxQ(int i,double x,double Q){
    double pdfs[14];
    xfxQarray(x, Q, pdfs);
    return pdfs[i+6];
  }

  void EvolutionHOPPET::xfxQarray(double x,double Q,double*pdfs){
    hoppetEval(x, Q, pdfs);
  }

  double EvolutionHOPPET::getAlphaS(double Q){
    return hoppetAlphaS(Q);
  }


  vector<double> EvolutionHOPPET::getXgrid() {
    vector<double> qx;
    qx.resize(2);
    qx[0] = exp(-1 * _ymax);
    qx[1] = 1.0;
    return qx;
  }

  vector<double> EvolutionHOPPET::getQgrid() {
    vector<double> qGrid;
    qGrid.resize(2);
    qGrid[0] = _Qmin;
    qGrid[1] = _Qmax;
    return qGrid;
  }
  
}

