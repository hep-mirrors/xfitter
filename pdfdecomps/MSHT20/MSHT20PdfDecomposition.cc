
#include "MSHT20PdfDecomposition.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include <iostream>
#include <iomanip>
#include <cmath>

namespace xfitter {

  /// the class factories, for dynamic loading
  extern "C" MSHT20PdfDecomposition* create(const char*name) {
    return new MSHT20PdfDecomposition(name);
  }

  // Constructor
  MSHT20PdfDecomposition::MSHT20PdfDecomposition(const char* inName) : BasePdfDecomposition(inName) {
  }

  const char*MSHT20PdfDecomposition::getClassName()const{return"MSHT20";}

  // Init at start:
  void MSHT20PdfDecomposition::atStart() {
    const YAML::Node node=XFITTER_PARS::getDecompositionNode(_name);
    //TODO: handle errors
    par_xuv  =getParameterisation(node["xuv"].as<string>());
    par_xdv  =getParameterisation(node["xdv"].as<string>());
    par_xs   =getParameterisation(node["xs"].as<string>());
    par_xr   =getParameterisation(node["xr"].as<string>());
    par_xsp  =getParameterisation(node["xsp"].as<string>());
    par_xsm  =getParameterisation(node["xsm"].as<string>());
    par_xg   =getParameterisation(node["xg"].as<string>());
    par_xgn  =getParameterisation(node["xgn"].as<string>());
  }

  void MSHT20PdfDecomposition::atIteration() {
    //Enforce sum rules
    // counting sum-rules for uv and dv
    par_xuv->setMoment(-1,2.0);
    par_xdv->setMoment(-1,1.0);
    // momentum sum-rule
    // quark part
    double xsumq=0;
    xsumq+=  par_xuv  ->moment(0);
    xsumq+=  par_xdv  ->moment(0);
    xsumq+=  par_xs   ->moment(0);
    // gluon part
    xsumq+=  par_xgn  ->moment(0);
    par_xg->setMoment(0,1-xsumq);
    // s-sbar=0 constraint
    // currently assume parametrization
    // s(x)-sbar(x)=A*(1-x)^eta*(1-x/x0)^delta
    // then int_0^1(s-sbar)=...
    // x0=delta/(delta+eta+1) [chatgpt]
    // TODO check
    // the order of parameters must be [A, eta, delta, x0]
    /*printf("dupa\n");fflush(stdout);
    double eta = par_xsm->getParameter(1);
    printf("dupa2\n");fflush(stdout);
    double delta = par_xsm->getParameter(2);
    double x0 = delta / (delta + eta + 1);
    printf("par_xsm->getNPar() = %d\n", par_xsm->getNPar());fflush(stdout);
    par_xsm->setParameter(par_xsm->getNPar()-1, x0);*/
  }

  // Returns a LHAPDF-style function, that returns PDFs in a physical basis for given x
  std::map<int,double>MSHT20PdfDecomposition::xfxMap(double x)const
  {
    double uv=(*par_xuv)(x);
    double dv=(*par_xdv)(x);
    double sp=(*par_xsp)(x);
    double sm=(*par_xsm)(x);
    double s=(sp+sm)/2.;
    double sbar=(sp-sm)/2.;
    double sea=(*par_xs)(x);
    //u+d=upd
    //d/u=dru -> d=u*dru
    //u+u*dru=upd -> u=upd/(1+dru)
    double ubar_plus_dbar=(sea-s-sbar)/2.;
    double dbar_over_ubar=(*par_xr)(x);
    double ubar=ubar_plus_dbar/(dbar_over_ubar+1);
    double dbar=ubar_plus_dbar-ubar;
    double g=(*par_xg)(x);
    double gn=(*par_xgn)(x);
    return{
      {-6,0},
      {-5,0},
      {-4,0},
      {-3,sbar},
      {-2,ubar},
      {-1,dbar},
      { 1,dbar+dv},
      { 2,ubar+uv},
      { 3,s},
      { 4,0},
      { 5,0},
      { 6,0},
      {21,g+gn}
    };
  }
}
