#include "Pion_FF_PdfParam.h"
#include "xfitter_cpp_base.h"
#include <cmath>
#include <cassert>
#include <iostream>
using namespace std;


namespace xfitter{
//for dynamic loading
extern"C" Pion_FF_PdfParam*create(const char*name){
  return new Pion_FF_PdfParam(name);
}
void Pion_FF_PdfParam::atStart(){
  using namespace std;
  BasePdfParam::atStart();
  const size_t n=getNPar();
  if(n<3){
    cerr<<"[ERROR] Too few parameters given to parameterisation \""<<_name<<"\", expected at least 3, got "<<n<<endl;
    hf_errlog(18120700,"F: Wrong number of parameters for a parameterisation, see stderr");
  }
  atIteration();
}
void Pion_FF_PdfParam::atIteration(){

}

double Pion_FF_PdfParam::operator()(double x)const{
  const unsigned int npar=getNPar();
  const double b=*pars[1]+2;
  const double c=*pars[2]+1;
  const double d=*pars[3];
  const double e=*pars[2] + *pars[4]+1;
  double norm = *pars[0]/(exp(lgamma(b)+lgamma(c)-lgamma(b+c)) + d*exp(lgamma(b)+lgamma(e)-lgamma(b+e)));

  double power= x*pow(x,(*pars[1]))*pow((1-x),(*pars[2]));
  double poly = 1;
  double xx = 1;
  poly= 1 + (*pars[3])*pow((1-x) , (*pars[4])) ;

  return norm*power*poly;
}

}
