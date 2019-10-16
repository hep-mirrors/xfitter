#include"FactorPdfParam.h"
#include"xfitter_cpp_base.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include<iostream>

namespace xfitter{
extern"C" FactorPdfParam*create(const char*s){return new FactorPdfParam(s);}

void FactorPdfParam::atStart(){
  using namespace std;
  YAML::Node node=XFITTER_PARS::getParameterisationNode(_name);

  //Get input
  try{
    input = xfitter::getParameterisation(node["input"].as<string>());
  }catch(YAML::TypedBadConversion<string>){
    cerr<<"[ERROR] Bad option \"input\" for parameterisation \""<<_name<<"\", expected name of another parameterisation"<<endl;
    hf_errlog(19070400, "F: Bad input for a Factor parameterisation, see stderr");
  }

  //Get factor
  try{
    factor = XFITTER_PARS::getParamD(node["factor"].as<string>());
  }catch(YAML::TypedBadConversion<string>){
    cerr<<"[ERROR] Bad option \"factor\" for parameterisation \""<<_name<<"\", expected parameter name"<<endl;
    hf_errlog(19070401, "F: Bad factor name for a Factor parameterisation, see stderr");
  }
}

double FactorPdfParam::operator()(double x)const{
  return (*factor) * (*input)(x);
}

double FactorPdfParam::moment(int n)const{
  return (*factor) * input->moment(n);
}

void FactorPdfParam::setMoment(int n,double val){
  input->setMoment(n, val / *factor );
}
}
