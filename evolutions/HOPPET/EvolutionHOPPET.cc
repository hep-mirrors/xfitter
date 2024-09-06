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
    
    int PtOrder = OrderMap(XFITTER_PARS::getParamS("Order"));
    // temporary: allow different orders in evolution and DIS SFs
    if (XFITTER_PARS::gParametersS.find("Order_HOPPET_Evolution") != XFITTER_PARS::gParametersS.end()) {
      PtOrder = OrderMap(XFITTER_PARS::getParamS("Order_HOPPET_Evolution"));
    }

    double dy = yamlNode["dy"].as<double>();
    hoppetStart(dy, PtOrder);
    int isFFNS = 0; // VFNS by default
    if(XFITTER_PARS::gParametersI.find("isFFNS") != XFITTER_PARS::gParametersI.end())
      isFFNS = XFITTER_PARS::gParametersI.at("isFFNS");
    if(isFFNS == 1) {
      int nflavour = XFITTER_PARS::gParametersI.at("NFlavour");
      hoppetSetFFN(nflavour);
    }
    else if(isFFNS == 0) {
      // TODO: check what will happen if these are free parameters (atConfigurationChange)
      const double* MCharm   = XFITTER_PARS::getParamD("mch");
      const double* MBottom  = XFITTER_PARS::getParamD("mbt");
      const double* MTop     = XFITTER_PARS::getParamD("mtp");
      hoppetSetVFN(*MCharm, *MBottom, *MTop);
    }
    else
    {
      hf_errlog(20240903, "F: Unsupported isFFNS = " + std::to_string(isFFNS));
    }
    atConfigurationChange();
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
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order"));
    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    const double* Q_ref      = XFITTER_PARS::getParamD("Mz");
    const double* Alphas_ref = XFITTER_PARS::getParamD("alphas");
    //const YAML::Node QGrid   = yamlNode["QGrid"];
    gPdfDecomp = _inPDFs;

    hoppetEvolve( *Alphas_ref, *Q_ref, PtOrder, 1.0, heralhc_init, *Q0);
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


  // Optional (can be done later):
  vector<double> EvolutionHOPPET::getXgrid() {
    hf_errlog(2024090401, "F: HOPPET getXgrid is not implemented yet");
    return vector<double>();
  }

  vector<double> EvolutionHOPPET::getQgrid() {
    hf_errlog(2024090402, "F: HOPPET getQgrid is not implemented yet");
    return vector<double>();
  }
  
}

