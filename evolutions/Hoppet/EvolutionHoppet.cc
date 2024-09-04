#include "EvolutionHoppet.h"
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
  extern "C" EvolutionHoppet*create(const char*name){
    return new EvolutionHoppet(name);
  }
  const char*EvolutionHoppet::getClassName()const{return "Hoppet";}

 void  heralhc_init(const double & x,
                    const double & Q,
                    double * pdf)                  
 {
/*
  double uv, dv;
  double ubar, dbar;
  double N_g=1.7, N_ls=0.387975;
  double N_uv=5.107200, N_dv=3.064320;
  double N_db=N_ls/2;

  uv = N_uv * pow(x,0.8) * pow((1-x),3);
  dv = N_dv * pow(x,0.8) * pow((1-x),4);
  dbar = N_db * pow(x,-0.1) * pow(1-x,6);
  ubar = dbar * (1-x);

  pdf[ 0+6] = N_g * pow(x,-0.1) * pow(1-x,5);
  pdf[-3+6] = 0.2*(dbar + ubar);
  pdf[ 3+6] = pdf[-3+6];
  pdf[ 2+6] = uv + ubar;
  pdf[-2+6] = ubar;
  pdf[ 1+6] = dv + dbar;
  pdf[-1+6] = dbar;

  pdf[ 4+6] = 0;
  pdf[ 5+6] = 0;
  pdf[ 6+6] = 0;
  pdf[-4+6] = 0;
  pdf[-5+6] = 0;
  pdf[-6+6] = 0;
*/	
  const std::map<int, int> ip =
      {{-5, -5},{-4, -4},{-3, -3}, {-2, -2}, {-1, -1}, {0, 21}, {1, 1}, {2, 2}, {3, 3},{4, 4},{5, 5}};
  for (auto &i : ip)
  {
    pdf[i.first + 6] = gPdfDecomp->xfxMap(x)[i.second]; 
    //double t = gPdfDecomp->xfxMap(x)[i.second]; 
    //std::cout << " from decomposition: " << i.first << " " << i.second << "  " << t << std::endl;
  }
 } 

void lha_unpolarized_dummy_pdf(const double &x, const double &Q, double *pdf)
{
    double uv, dv;
    double ubar, dbar;
    double N_g = 1.7, N_ls = 0.387975;
    double N_uv = 5.107200, N_dv = 3.064320;
    double N_db = N_ls / 2;

    uv = N_uv * pow(x, 0.8) * pow((1 - x), 3);
    dv = N_dv * pow(x, 0.8) * pow((1 - x), 4);
    dbar = N_db * pow(x, -0.1) * pow(1 - x, 6);
    ubar = dbar * (1 - x);

    pdf[0 + 6] = N_g * pow(x, -0.1) * pow(1 - x, 5);
    pdf[-3 + 6] = 0.2 * (dbar + ubar);
    pdf[3 + 6] = pdf[-3 + 6];
    pdf[2 + 6] = uv + ubar;
    pdf[-2 + 6] = ubar;
    pdf[1 + 6] = dv + dbar;
    pdf[-1 + 6] = dbar;

    pdf[4 + 6] = 0;
    pdf[5 + 6] = 0;
    pdf[6 + 6] = 0;
    pdf[-4 + 6] = 0;
    pdf[-5 + 6] = 0;
    pdf[-6 + 6] = 0;
}

  void EvolutionHoppet::atStart()
  {
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    _inPDFs=XFITTER_PARS::getInputDecomposition(yamlNode);
    const YAML::Node xGrid = yamlNode["xGrid"];
    
    //const int PtOrder = OrderMap(XFITTER_PARS::getParamS("Order")) - 1; // here was -1
    const int PtOrder = OrderMap(XFITTER_PARS::getParamS("Order")); // here was -1
    double dy = 0.1;
    hoppetStart(dy, PtOrder);
    int isFFNS = 0; // VFNS by default
    if(XFITTER_PARS::gParametersI.find("isFFNS") != XFITTER_PARS::gParametersI.end())
      isFFNS = XFITTER_PARS::gParametersI.at("isFFNS");
    if(isFFNS == 1) {
      int nflavour = XFITTER_PARS::gParametersI.at("NFlavour");
      hoppetSetFFN(nflavour);
    }
    else if(isFFNS == 0) {
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

  
    void EvolutionHoppet::atIteration()
  {
    // use  https://github.com/hoppet-code/hoppet/blob/master/example_f77/cpp_tabulation_example.cc

    
    std::cout << " HERE WE ARE in HOPPET " << std::endl;
    //exit(0);
    
    const YAML::Node yamlNode=XFITTER_PARS::getEvolutionNode(_name);
    // Retrieve the relevant parameters needed to compute the evolutions
    //const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order")) - 1;//here was -1
    const int     PtOrder    = OrderMap(XFITTER_PARS::getParamS("Order"));//here was -1
    const double* Q0         = XFITTER_PARS::getParamD("Q0");
    const double* Q_ref      = XFITTER_PARS::getParamD("Mz");
    const double* Alphas_ref = XFITTER_PARS::getParamD("alphas");
    const YAML::Node QGrid   = yamlNode["QGrid"];
    gPdfDecomp = _inPDFs;
    // add hoppetEvolve
     //double asQ0 = 0.118;
    const bool param_coefs = true;
    const double xmuR = 1.;
    const double xmuF = 1.;

    double Qmax = 13000.0;
    double Qmin = 1.0;
    int order = -6;
    double ymax = 16.0;
    double dy = 0.05;
    double dlnlnQ = dy / 4.0;
    int nloop = 3;
    double minQval = min(xmuF * Qmin, Qmin);
    double maxQval = max(xmuF * Qmax, Qmax);
  
   //hoppetStartExtended(ymax, dy, minQval, maxQval, dlnlnQ, nloop, order, factscheme_MSbar);
    
    int nflav = -5;
    int order_max = 4;
    int sc_choice = scale_choice_Q;
    double zmass = 91.1876;
    double wmass = 80.377;
    
   // hoppetStartStrFctExtended(order_max, nflav, scale_choice_Q, zmass, param_coefs, wmass, zmass);
   
    hoppetEvolve( *Alphas_ref, *Q_ref, PtOrder, 1.0, heralhc_init, *Q0);
    double f[13];
    hoppetEval(0.001,10.,f);
    printf("f[6] = %f  as = %f\n", f[6], hoppetAlphaS(10.));
    std::cout << " HERE WE ARE OUT OF HOPPET " << std::endl;
    
    
   
    
    
    
    double asQ = 0.35;
    double QH = sqrt(2.0);
    double muR_Q = 1.0;
      
    // Evolve the PDF
   //hoppetEvolve(asQ, QH, nloop, muR_Q, lha_unpolarized_dummy_pdf, QH);

    // Initialize structure functions(we dk need it or not)
   // hoppetInitStrFct(order_max, param_coefs, xmuR, xmuF);
    
  }
  
  
  std::map<int,double>EvolutionHoppet::xfxQmap(double x,double Q){
    std::cout << " HERE WE ARE in HOPPET A " << std::endl;
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

  double EvolutionHoppet::xfxQ(int i,double x,double Q){
    double pdfs[14];
    xfxQarray(x, Q, pdfs);
    return pdfs[i+6];
  }

  void EvolutionHoppet::xfxQarray(double x,double Q,double*pdfs){
    hoppetEval(x, Q, pdfs);
  }

  double EvolutionHoppet::getAlphaS(double Q){
    return hoppetAlphaS(Q);
  }


  // Optional (can be done later):
  vector<double> EvolutionHoppet::getXgrid() {
  }

  vector<double> EvolutionHoppet::getQgrid() {
  }
  
 
}

