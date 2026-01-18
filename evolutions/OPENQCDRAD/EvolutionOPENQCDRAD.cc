#include "EvolutionOPENQCDRAD.h"
#include "xfitter_pars.h"
#include "xfitter_cpp_base.h"
#include "ABM.h"
//#include <iostream>
//#include <cmath>

using namespace std;

// Global var to hold current pdfDecomposition
xfitter::BasePdfDecomposition *gPdfDecomp = nullptr;
extern "C" {
  double xqg0_(const int& k, const int& iq, const double& xb, const int& ix) {
    std::map<int, int> ip =
      {{1, 21},{2, 1},{3, -1}, {4, 2}, {5, -2}, {6, 3}, {7, -3}};
    if(ip.find(iq) != ip.end()) {
      return gPdfDecomp->xfxMap(xb)[ip[iq]];
    }
    else {
      throw 42;
      return 0;
    }
  };
}

namespace xfitter
{

  // the class factories
  extern "C" EvolutionOPENQCDRAD*create(const char*name){
    return new EvolutionOPENQCDRAD(name);
  }
  const char*EvolutionOPENQCDRAD::getClassName()const{return "OPENQCDRAD";}

  void EvolutionOPENQCDRAD::atStart()
  {
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    gPdfDecomp=XFITTER_PARS::getInputDecomposition(yamlNode);
    //const YAML::Node xGrid = yamlNode["xGrid"];
    
    _order = OrderMap(XFITTER_PARS::getParamS("Order"));
    // temporary: allow different orders in evolution and DIS SFs
    if (XFITTER_PARS::gParametersS.find("Order_OPENQCDRAD_Evolution") != XFITTER_PARS::gParametersS.end()) {
      _order = OrderMap(XFITTER_PARS::getParamS("Order_OPENQCDRAD_Evolution"));
    }

    //double dy = yamlNode["dy"].as<double>();
    //OPENQCDRADStart(dy, _order);
    //double dy_over_dlnlnQ = yamlNode["dy_over_dlnlnQ"].as<double>();
    //double dlnlnQ = dy/dy_over_dlnlnQ;
    //_ymax = yamlNode["ymax"].as<double>();
    //_Qmin = yamlNode["Qmin"].as<double>();
    //_Qmax = yamlNode["Qmax"].as<double>();
    //int order_interpol = yamlNode["order_interpol"].as<int>();
    //OPENQCDRADStartExtended(_ymax, dy, _Qmin, _Qmax, dlnlnQ, _order, order_interpol, factscheme_MSbar);
    //_isFFNS = 0; // VFNS by default
    //if (yamlNode["isFFNS"]) {
    //  _isFFNS = yamlNode["isFFNS"].as<int>();
    //}
    //else {
    //  if(XFITTER_PARS::gParametersI.find("isFFNS") != XFITTER_PARS::gParametersI.end()) {
    //    _isFFNS = XFITTER_PARS::gParametersI.at("isFFNS");
    //  }
    //}
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
    //if(_isFFNS == 1) {
      //OPENQCDRADSetFFN(_nflavour);
    //}
    //else if(_isFFNS == 0) {
      _mch = XFITTER_PARS::getParamD("mch");
      _mbt = XFITTER_PARS::getParamD("mbt");
      //_mtp = XFITTER_PARS::getParamD("mtp");
      //OPENQCDRADSetVFN(*_mch, *_mbt, *_mtp);
      //OPENQCDRADSetPoleMassVFN(*_mch, *_mbt, *_mtp);
      //if (_msbar) {
      //  OPENQCDRADSetMSbarMassVFN(*_mch, *_mbt, *_mtp);
      //}  }
    //}
    //else
    //{
    //  hf_errlog(20240903, "F: Unsupported isFFNS = " + std::to_string(_isFFNS));
    //}
    atConfigurationChange();
  }
  
    void EvolutionOPENQCDRAD::atIteration()
  {    
    /*const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
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
    gPdfDecomp = _inPDFs;*/

    // based on pdfs/pdfevol.f from OPENQCDRAD benchmark
    abm::set_scheme_and_order_for_evolution(*_mch, *_mbt);
  }
  
  
  std::map<int,double>EvolutionOPENQCDRAD::xfxQmap(double x,double Q){
    double pdfs[14] = {0.0};
    xfxQarray(x, Q, pdfs);
    std::map<int, double> res;
    
    const int npdfMax = 6;

    for (int ipdf = -6; ipdf <= npdfMax; ipdf++)
      {
	int ii = (ipdf == 0) ? 21 : ipdf;
	// photon PDF (not implemented):
	if (ipdf == 7)
	  ii = 22;
	res[ii] = pdfs[ipdf+6];
      }
    return res;
  }

  double EvolutionOPENQCDRAD::xfxQ(int i,double x,double Q){
    double pdfs[14] = {0.0};
    //printf("xfxQ(x=%.1e) = ", x); for (int i = 0; i < 14; i++) printf(" [%d]%+.1e", i, pdfs[i]); printf("\n");
    xfxQarray(x, Q, pdfs);
    return pdfs[i+6];
  }

  void EvolutionOPENQCDRAD::xfxQarray(double x,double Q,double*pdfs){
    const int kp = 0;
    for (int i = 0; i < 13; i++) {
      pdfs[i] = 0.;
      //pdfs[i] = abm::xqg(i, x, Q, kp);
      //pdfs[i] = abm::xqg(i+1, x, Q, kp);
      //pdfs[i] = abm::xqg(kp, x, Q, i);
    }
    const double Q2 = Q*Q;
    pdfs[6] = abm::xqg(1, x, Q2, kp);
    pdfs[5] = abm::xqg(3, x, Q2, kp);
    pdfs[4] = abm::xqg(5, x, Q2, kp);
    pdfs[3] = abm::xqg(7, x, Q2, kp);
    pdfs[7] = abm::xqg(2, x, Q2, kp);
    pdfs[8] = abm::xqg(4, x, Q2, kp);
    pdfs[9] = abm::xqg(6, x, Q2, kp);
    //pdfs[7] = pdfs[5] + abm::xqg(2, x, Q2, kp);
    //pdfs[8] = pdfs[4] + abm::xqg(4, x, Q2, kp);
    //pdfs[9] = pdfs[3] + abm::xqg(6, x, Q2, kp);
  }

  double EvolutionOPENQCDRAD::getAlphaS(double Q){
    const double Q2 = Q*Q;
    return abm::alphas(Q2);
  }

  /*vector<double> EvolutionOPENQCDRAD::getXgrid() {
    vector<double> qx;
    qx.resize(2);
    qx[0] = exp(-1 * _ymax);
    qx[1] = 1.0;
    return qx;
  }

  vector<double> EvolutionOPENQCDRAD::getQgrid() {
    vector<double> qGrid;
    qGrid.resize(2);
    qGrid[0] = _Qmin;
    qGrid[1] = _Qmax;
    return qGrigetXgrid;
  }*/
  
}

