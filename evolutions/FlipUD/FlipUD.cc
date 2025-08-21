#include"FlipUD.h"
#include"xfitter_cpp_base.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"iostream"
#include<utility> //for swap
using namespace std;
namespace xfitter{
extern"C" FlipUD_Evol*create(const char*s){return new FlipUD_Evol(s);}//for dynamic loading
void FlipUD_Evol::atStart(){atConfigurationChange();}
void FlipUD_Evol::atConfigurationChange(){
  const YAML::Node node=XFITTER_PARS::getEvolutionNode(_name)["input"];
  try{
    input=get_evolution(node.as<string>());
  }catch(const YAML::BadConversion&ex){
    if(!node){
      //If no input is provided, use default evolution
      input=defaultEvolutionInstance();
      cerr<<"[INFO] Evolution \""<<_name<<"\" of class "<<getClassName()<<": \"input\" not provided; using default evolution \""<<input->_name<<"\""<<endl;
    }else{
      cerr<<"[ERROR] In evolution "<<_name<<" of class "<<getClassName()<<": failed to convert input to string, YAML node:"<<node<<endl;
      hf_errlog(18122501,"F: Bad input to evolution FlipCharge, see stderr");
    }
  }
  _flip = XFITTER_PARS::getEvolutionNode(_name)["flip"] ? XFITTER_PARS::getEvolutionNode(_name)["flip"].as<double>() : 1.0;
}
map<int,double>FlipUD_Evol::xfxQmap(double x,double Q){
  double p[13];
  input->xfxQarray(x,Q,p);
  return{
    {-6,p[ 0]},
    {-5,p[ 1]},
    {-4,p[ 2]},
    {-3,p[ 3]},
    {-2,p[ 5]*_flip+p[4]*(1-_flip)},
    {-1,p[ 4]*_flip+p[5]*(1-_flip)},
    {21,p[ 6]},
    { 1,p[ 8]*_flip+p[7]*(1-_flip)},
    { 2,p[ 7]*_flip+p[8]*(1-_flip)},
    { 3,p[ 9]},
    { 4,p[10]},
    { 5,p[11]},
    { 6,p[12]}
  };
}
void FlipUD_Evol::xfxQarray(double x,double Q,double*p){
  input->xfxQarray(x,Q,p);
  //swap(p[4],p[5]);
  //swap(p[7],p[8]);
  double p4 = p[4];
  p[4] = p[5]*_flip+p4*(1-_flip);
  p[5] = p4*_flip+p[5]*(1-_flip);
  double p7 = p[7];
  p[7] = p[8]*_flip+p7*(1-_flip);
  p[8] = p7*_flip+p[8]*(1-_flip);
}
double FlipUD_Evol::xfxQ(int i,double x,double Q){
  switch(i){
    case  1:return input->xfxQ( 2,x,Q)*_flip+input->xfxQ(1,x,Q)*(1-_flip);
    case -1:return input->xfxQ(-2,x,Q)*_flip+input->xfxQ(-1,x,Q)*(1-_flip);
    case  2:return input->xfxQ( 1,x,Q)*_flip+input->xfxQ(2,x,Q)*(1-_flip);
    case -2:return input->xfxQ(-1,x,Q)*_flip+input->xfxQ(-2,x,Q)*(1-_flip);
    default:return input->xfxQ( i,x,Q);
  }
}
double FlipUD_Evol::getAlphaS(double Q){return input->getAlphaS(Q);}

vector<double> FlipUD_Evol::getXgrid(){
  return input->getXgrid();
}

vector<double> FlipUD_Evol::getQgrid(){
  return input->getQgrid();
}
}
