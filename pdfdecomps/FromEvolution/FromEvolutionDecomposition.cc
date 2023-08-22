 
/*
   @file FromEvolutionDecomposition.cc
   @date 2023-08-18
   @author  Sasha Zenaiev [oleksandr.zenaiev@desy.de]
*/

#include"FromEvolutionDecomposition.h"
#include"xfitter_pars.h"
#include"xfitter_cpp_base.h"
#include"xfitter_steer.h"
#include <iostream>

namespace xfitter {
  //For dynamic loading:
  extern "C" FromEvolutionDecomposition*create(const char*name){
    return new FromEvolutionDecomposition(name);
  }
  //_________________________________________________________________________________
  FromEvolutionDecomposition::FromEvolutionDecomposition(const char*name):BasePdfDecomposition{name}{}
  FromEvolutionDecomposition::~FromEvolutionDecomposition(){}
  const char*FromEvolutionDecomposition::getClassName()const{return"FromEvolution";}
  void FromEvolutionDecomposition::atStart(){
    YAML::Node pars=XFITTER_PARS::getDecompositionNode(_name);
    try{
      _Q0 = *XFITTER_PARS::getParamD(pars["Q0"] ? pars["Q0"].as<string>() : "Q0");
    }catch(YAML::TypedBadConversion<std::string>&ex){
      std::cerr<<"[ERROR] Bad Q0 parameter name given for FromEvolution decomposition "<<_name<<":\n"<<pars<<"\n[/ERROR]"<<std::endl;
      hf_errlog(18090310,"F: In FromEvolutionDecomposition: failed to convert YAML node \"Q0\" to string, see stderr");
    }
  }
  void FromEvolutionDecomposition::atIteration(){
    // need to read evolution name at iteration rather than at start, because at start evolutions do not exist yet
    YAML::Node pars=XFITTER_PARS::getDecompositionNode(_name);
    try{
      _evolution=get_evolution((pars["evolution"].as<std::string>()));
    }catch(YAML::TypedBadConversion<std::string>&ex){
      std::cerr<<"[ERROR] Bad evolution name given for FromEvolution decomposition "<<_name<<":\n"<<pars<<"\n[/ERROR]"<<std::endl;
      hf_errlog(18090310,"F: In FromEvolutionDecomposition: failed to convert YAML node \"evolution\" to string, see stderr");
    }
  }
  std::map<int,double>FromEvolutionDecomposition::xfxMap(double x)const{
    return _evolution->xfxQmap(x,  _Q0);
  }
}



