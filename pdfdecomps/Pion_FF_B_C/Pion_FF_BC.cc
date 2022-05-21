 
/*
   @file Pion-FF_PdfDecomposition.cc
   @date 2018-08-14
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-11
*/

#include "Pion_FF_BC.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;


namespace xfitter
{
  /// the class factories
  extern "C" Pion_FF_BC*create(const char*name){
    return new Pion_FF_BC(name);
  }

  
  //_________________________________________________________________________________
  Pion_FF_BC::Pion_FF_BC(const char*name):BasePdfDecomposition{name}{}
  const char*Pion_FF_BC::getClassName()const{return"Pion_FF_BC";}

  //_________________________________________________________________________________
  void Pion_FF_BC::atStart(){
    const YAML::Node node=XFITTER_PARS::getDecompositionNode(_name);
    //TODO: handle errors
    par_xup  =getParameterisation(node["xup"].as<string>());
  //  par_xdp  =getParameterisation(node["xdp"].as<string>());
    par_xcp  =getParameterisation(node["xcp"].as<string>());
    par_xbp =getParameterisation(node["xbp"].as<string>());
    par_xsp   =getParameterisation(node["xsp"].as<string>());
    par_xg   =getParameterisation(node["xg"].as<string>());
  }

  void Pion_FF_BC::atIteration() {
    //Enforce sum rules
    // counting sum-rules for uv and dv
   // par_xuv->setMoment(-1,2.0);
   // par_xdv->setMoment(-1,1.0);
    // momentum sum-rule
    // quark part
  //  double xsumq=0;
   // xsumq+=  par_xuv  ->moment(0);
  //  xsumq+=  par_xdv  ->moment(0);
///xsumq+=2*par_xubar->moment(0);
  //  xsumq+=2*par_xdbar->moment(0);
  //  xsumq+=2*par_xs   ->moment(0);
    // gluon part
 //   par_xg->setMoment(0,1-xsumq);
  // cout << "Hello HamedpionFF!";
  }
  std::map<int,double>Pion_FF_BC::xfxMap(double x)const
  {
    double bbar=0.5*(*par_xbp)(x);
    double cbar=0.5*(*par_xcp)(x);
    double ubar=0.5*(*par_xup)(x);
   // double dbar=(*par_xdbar)(x);
  //  double u=(*par_xuv)(x)+ubar;
  //  double d=(*par_xdv)(x)+dbar;
    double s=0.5*(*par_xsp)(x);
    double g=(*par_xg)(x);
    return{
      {-6,0},
      {-5,bbar},
      {-4,cbar},
      {-3,s},
      {-2,ubar},
      {-1,ubar},
      { 1,ubar},
      { 2,ubar},
      { 3,s},
      { 4,cbar},
      { 5,bbar},
      { 6,0},
      {21,g}
    };
  }
}


