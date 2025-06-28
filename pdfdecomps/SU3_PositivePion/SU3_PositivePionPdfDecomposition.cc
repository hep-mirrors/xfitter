/*
   @file SU3_PositivePionPdfDecomposition.cc
   @date 2018-08-07
   @author  AddPdfDecomposition.py
   Created by  AddPdfDecomposition.py on SU3_PositivePion
*/
#include"xfitter_cpp_base.h"
#include<iostream>
#include"SU3_PositivePionPdfDecomposition.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
using namespace std;
//For dynamic loading:
namespace xfitter{
extern "C" SU3_PositivePionPdfDecomposition*create(const char*name){
  return new SU3_PositivePionPdfDecomposition(name);
}
SU3_PositivePionPdfDecomposition::SU3_PositivePionPdfDecomposition(const char*name):BasePdfDecomposition{name}{}
const char*SU3_PositivePionPdfDecomposition::getClassName()const{return"SU3_PositivePion";}
BasePdfParam*getParam(const BasePdfDecomposition*self,const YAML::Node&node,const char*s){
  try{
    return getParameterisation(node[s].as<string>());
  }catch(const YAML::InvalidNode&ex){
    if(node[s].IsNull()){
      cerr<<"[ERROR] No \""<<s<<"\" parameterisation given for decomposition \""<<self->_name<<"\""<<endl;
      hf_errlog(18092410,"F: Error in decomposition parameters, details written to stderr");
    }else throw ex;
  }catch(const YAML::BadConversion&ex){
    cerr<<"[ERROR] Bad parameter \""<<s<<"\" given for decomposition \""<<self->_name<<"\""<<endl;
    hf_errlog(18092410,"F: Error in decomposition parameters, details written to stderr");
  }
  return nullptr;//unreachable, suppress warning
}
// Init at start:
void SU3_PositivePionPdfDecomposition::atStart(){
  const YAML::Node node=XFITTER_PARS::getDecompositionNode(_name);
  //get parameterisation usually doesn't throw
  par_v=getParam(this,node,"valence");
  par_S=getParam(this,node,"sea");
  par_g=getParam(this,node,"gluon");
}
void SU3_PositivePionPdfDecomposition::atIteration() {
  //Enforce sum rules
  //Valence sum
  par_v->setMoment(-1,2);
  //Momentum sum
  par_g->setMoment(0,1-par_S->moment(0)-par_v->moment(0));
}
map<int,double>SU3_PositivePionPdfDecomposition::xfxMap(const double x)const{
  const double v=(*par_v)(x);
  const double S=(*par_S)(x);
  const double g=(*par_g)(x);
  const double d=S * (1./6.);
  const double u=v * 0.5 + d;
  return{
    {-6,0},
    {-5,0},
    {-4,0},
    {-3,d},//sbar
    {-2,d},//ubar
    {-1,u},//dbar
    { 1,d},//d
    { 2,u},//u
    { 3,d},//s
    { 4,0},
    { 5,0},
    { 6,0},
    {21,g}
  };
}
}
