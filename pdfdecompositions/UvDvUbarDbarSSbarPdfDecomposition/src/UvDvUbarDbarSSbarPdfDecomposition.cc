
/*
   @file UvDvUbarDbarSSbarPdfDecomposition.cc
   @date 2019-02-25
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2019-02-25
*/

#include "UvDvUbarDbarSSbarPdfDecomposition.h"
#include "xfitter_pars.h"
#include "xfitter_steer.h"
#include <iostream>
#include <iomanip>
#include <cmath>

namespace xfitter {

  /// the class factories, for dynamic loading
  extern "C" UvDvUbarDbarSSbarPdfDecomposition* create(const char*name) {
    return new UvDvUbarDbarSSbarPdfDecomposition(name);
  }

  // Constructor
  UvDvUbarDbarSSbarPdfDecomposition::UvDvUbarDbarSSbarPdfDecomposition(const char* inName) : BasePdfDecomposition(inName) {
  }

  const char*UvDvUbarDbarSSbarPdfDecomposition::getClassName()const{return"UvDvUbarDbarSSbar";}

  // Init at start:
  void UvDvUbarDbarSSbarPdfDecomposition::atStart() {
    const YAML::Node node=XFITTER_PARS::getDecompositionNode(_name);
    //TODO: handle errors
    par_xuv  =getParameterisation(node["xuv"].as<string>());
    par_xdv  =getParameterisation(node["xdv"].as<string>());
    par_xubar=getParameterisation(node["xubar"].as<string>());
    par_xdbar=getParameterisation(node["xdbar"].as<string>());
    par_xs   =getParameterisation(node["xs"].as<string>());
    par_xsbar=getParameterisation(node["xsbar"].as<string>());
    par_xg   =getParameterisation(node["xg"].as<string>());
  }

  void UvDvUbarDbarSSbarPdfDecomposition::atIteration() {
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
    xsumq+=  par_xs   ->moment(0);
    xsumq+=  par_xsbar->moment(0);
    // gluon part
    par_xg->setMoment(0,1-xsumq);
  }

  // Returns a LHAPDF-style function, that returns PDFs in a physical basis for given x
  std::map<int,double>UvDvUbarDbarSSbarPdfDecomposition::xfxMap(double x)const
  {
    double ubar=(*par_xubar)(x);
    double dbar=(*par_xdbar)(x);
    double u=(*par_xuv)(x)+ubar;
    double d=(*par_xdv)(x)+dbar;
    double s=(*par_xs)(x);
    double sbar=(*par_xsbar)(x);
    double g=(*par_xg)(x);
    return{
      {-6,0},
      {-5,0},
      {-4,0},
      {-3,sbar},
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
  }
}
