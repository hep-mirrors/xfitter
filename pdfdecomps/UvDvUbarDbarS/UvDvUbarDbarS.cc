 
/*
   @file UvDvUbarDbarS_PdfDecomposition.cc
   @date 2018-08-14
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-11
*/

#include "UvDvUbarDbarS.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include <iostream>
#include <iomanip>
#include <cmath>

namespace xfitter
{
  /// the class factories
  extern "C" UvDvUbarDbarS*create(const char*name){
    return new UvDvUbarDbarS(name);
  }

  
  //_________________________________________________________________________________
  UvDvUbarDbarS::UvDvUbarDbarS(const char*name):BasePdfDecomposition{name}{}
  const char*UvDvUbarDbarS::getClassName()const{return"UvDvUbarDbarS";}

  //_________________________________________________________________________________
  void UvDvUbarDbarS::atStart(){
    const YAML::Node node=XFITTER_PARS::getDecompositionNode(_name);
    //TODO: handle errors
    par_xuv  =getParameterisation(node["xuv"].as<string>());
    par_xdv  =getParameterisation(node["xdv"].as<string>());
    par_xubar=getParameterisation(node["xubar"].as<string>());
    par_xdbar=getParameterisation(node["xdbar"].as<string>());
    par_xs   =getParameterisation(node["xs"].as<string>());
    par_xg   =getParameterisation(node["xg"].as<string>());

    //try optinal photon:
    if ( node["xgamma"] ) 
      par_xgamma = getParameterisation(node["xgamma"].as<string>());
  }

  void UvDvUbarDbarS::atIteration() {
    //Enforce sum rules
    // counting sum-rules for uv and dv
    par_xuv->setMoment(-1,2.0);
    par_xdv->setMoment(-1,1.0);
    // momentum sum-rule
    // quark part
    double xsumq=0;
    xsumq+=  par_xuv  ->moment(0);
    xsumq+=  par_xdv  ->moment(0);
    xsumq+=2*par_xubar->moment(0);
    xsumq+=2*par_xdbar->moment(0);
    xsumq+=2*par_xs   ->moment(0);
    // maybe gamma
    if (par_xgamma != nullptr) 
      xsumq += par_xgamma -> moment(0);
    // gluon part
    par_xg->setMoment(0,1-xsumq);
  }
  std::map<int,double>UvDvUbarDbarS::xfxMap(double x)const
  {
    double ubar=(*par_xubar)(x);
    double dbar=(*par_xdbar)(x);
    double u=(*par_xuv)(x)+ubar;
    double d=(*par_xdv)(x)+dbar;
    double s=(*par_xs)(x);
    double g=(*par_xg)(x);
    std::map<int,double> out = {
      {-6,0},
      {-5,0},
      {-4,0},
      {-3,s},
      {-2,ubar},
      {-1,dbar},
      { 1,d},
      { 2,u},
      { 3,s},
      { 4,0},
      { 5,0},
      { 6,0},
      {21,g}
    };
    if (par_xgamma != nullptr) {
      out[22] = (*par_xgamma)(x);
    }
    return out;
  }
}


