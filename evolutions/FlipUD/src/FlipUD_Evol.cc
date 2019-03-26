#include"FlipUD_Evol.h"
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
  }catch(const YAML::InvalidNode&ex){
    if(node.IsNull()){
      cerr<<"[ERROR] In evolution "<<_name<<" of class "<<getClassName()<<": no input evolution"<<endl;
      hf_errlog(18122510,"F: No input to evolution FlipUD, see stderr");
    }else throw ex;
  }catch(const YAML::BadConversion&ex){
    cerr<<"[ERROR] In evolution "<<_name<<" of class "<<getClassName()<<": failed to convert input to string, YAML node:"<<node<<endl;
    hf_errlog(18122511,"F: Bad input to evolution FlipUD, see stderr");
  }
}
function<map<int,double>(double const&x,double const&Q)>FlipUD_Evol::xfxQMap(){
  auto f=input->xfxQMap();
  return [f](double const&x,double const&Q)->map<int,double>{
    map<int,double>M=f(x,Q);
    swap(M.at( 1),M.at( 2));
    swap(M.at(-1),M.at(-2));
    return M;
  };
}
function<void(double const&x,double const&Q,double*pdfs)>FlipUD_Evol::xfxQArray(){
  auto f=input->xfxQArray();
  return [f](double const&x,double const&Q,double*p){
    f(x,Q,p);
    swap(p[4],p[5]);
    swap(p[7],p[8]);
  };
}
function<double(int const&i,double const&x,double const&Q)>FlipUD_Evol::xfxQDouble(){
  auto f=input->xfxQDouble();
  return [f](const int&i,const double&x,const double&Q)->double{
    switch(i){
      case  1:return f( 2,x,Q);
      case -1:return f(-2,x,Q);
      case  2:return f( 1,x,Q);
      case -2:return f(-1,x,Q);
      default:return f(i,x,Q);
    }
  };
}
function<double(double const&Q)>FlipUD_Evol::AlphaQCD(){return input->AlphaQCD();}
}
