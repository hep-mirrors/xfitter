#include"FlipChargeEvol.h"
#include"xfitter_cpp_base.h"
#include"xfitter_pars.h"
#include"xfitter_steer.h"
#include"iostream"
#include<utility> //for swap
using namespace std;
namespace xfitter{
extern"C" FlipChargeEvol*create(const char*s){return new FlipChargeEvol(s);}//for dynamic loading
void FlipChargeEvol::atStart(){atConfigurationChange();}
void FlipChargeEvol::atConfigurationChange(){
  const YAML::Node node=XFITTER_PARS::getEvolutionNode(_name)["input"];
  try{
    input=get_evolution(node.as<string>());
  }catch(const YAML::InvalidNode&ex){
    if(node.IsNull()){
      cerr<<"[ERROR] In evolution "<<_name<<" of class "<<getClassName()<<": no input evolution"<<endl;
      hf_errlog(18122500,"F: No input to evolution FlipCharge, see stderr");
    }else throw ex;
  }catch(const YAML::BadConversion&ex){
    cerr<<"[ERROR] In evolution "<<_name<<" of class "<<getClassName()<<": failed to convert input to string, YAML node:"<<node<<endl;
    hf_errlog(18122501,"F: Bad input to evolution FlipCharge, see stderr");
  }
}
function<map<int,double>(double const&,double const&)>FlipChargeEvol::xfxQMap(){
  auto f=input->xfxQMap();
  return [f](double const&x,double const&Q)->map<int,double>{
    map<int,double>M=f(x,Q);//I don't know how to avoid copying a map with current interface --Ivan
    swap(M.at(1),M.at(-1));
    swap(M.at(2),M.at(-2));
    swap(M.at(3),M.at(-3));
    swap(M.at(4),M.at(-4));
    swap(M.at(5),M.at(-5));
    swap(M.at(6),M.at(-6));
    return M;
  };
}
function<void(double const&x,double const&Q,double*pdfs)>FlipChargeEvol::xfxQArray(){
  auto f=input->xfxQArray();
  return [f](double const&x,double const&Q,double*p){
    f(x,Q,p);
    swap(p[0],p[12]);
    swap(p[1],p[11]);
    swap(p[2],p[10]);
    swap(p[3],p[ 9]);
    swap(p[4],p[ 8]);
    swap(p[5],p[ 7]);
  };
}
function<double(int const&i,double const&x,double const&Q)>FlipChargeEvol::xfxQDouble(){
  auto f=input->xfxQDouble();
  return [f](const int&i,const double&x,const double&Q)->double{return f(-i,x,Q);};
}
function<double(double const&Q)>FlipChargeEvol::AlphaQCD(){return input->AlphaQCD();}
}
