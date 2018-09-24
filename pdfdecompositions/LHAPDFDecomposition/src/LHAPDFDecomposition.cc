 
/*
   @file LHAPDF.cc
   @date 2018-07-12
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on 2018-07-12
*/

#include"LHAPDFDecomposition.h"
#include"xfitter_pars.h"
#include"xfitter_cpp_base.h"
#include"CheckForPDF.h"

namespace xfitter {
  //For dynamic loading:
  extern "C" LHAPDFDecomposition*create(const char*name){
    return new LHAPDFDecomposition(name);
  }
  //_________________________________________________________________________________
  LHAPDFDecomposition::LHAPDFDecomposition(const char*name):BasePdfDecomposition{name}{}
  LHAPDFDecomposition::~LHAPDFDecomposition(){if(_pdf)delete _pdf;}
  const char*LHAPDFDecomposition::getClassName()const{return"LHAPDF";}
  void LHAPDFDecomposition::initAtStart(){
    YAML::Node pars=XFITTER_PARS::getDecompositionNode(_name);
    string setName;
    int member;
    try{
      setName=pars["set"].as<std::string>();
    }catch(YAML::TypedBadConversion<std::string>&ex){
      hf_errlog(18090310,"F: In LHAPDFDecomposition::initAtStart: failed to convert YAML node \"set\" to string; printing node to stderr");
      std::cerr<<pars<<std::endl;
    }

    // check if exists first
    CheckForPDF(setName.c_str());
    member=pars["member"].as<int>();
    if(_pdf)delete _pdf;
    _pdf=LHAPDF::mkPDF(setName,member);
  }

  //_________________________________________________________________________________
  std::function<std::map<int,double>(const double& x)> LHAPDFDecomposition::f0() const
  {
    return [=] (const double& x)->std::map<int,double>
      {
        return _pdf->xfxQ(x, *(XFITTER_PARS::gParameters.at("Q0")));
      };
  }
}



