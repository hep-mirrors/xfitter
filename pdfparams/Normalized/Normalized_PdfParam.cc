#include "Normalized_PdfParam.h"
#include "xfitter_cpp_base.h"
#include <cmath>
#include <cassert>
#include <iostream>

namespace xfitter{
//for dynamic loading
extern"C" Normalized_PdfParam*create(const char*name){
  return new Normalized_PdfParam(name);
}
void Normalized_PdfParam::atStart(){
  using namespace std;
  BasePdfParam::atStart();
  const size_t n=getNPar();
  if(n<3){
    cerr<<"[ERROR] Too few parameters given to parameterisation \""<<_name<<"\", expected at least 3, got "<<n<<endl;
    hf_errlog(18120700,"F: Wrong number of parameters for a parameterisation, see stderr");
  }
  atIteration();
}
void Normalized_PdfParam::atIteration(){
  updateNormalization();
}
void Normalized_PdfParam::updateNormalization(){
  //norm=A/beta(B+1,C+1)
  const double b=*pars[1]+1;
  const double c=*pars[2]+1;
  norm=(*pars[0])*exp(-lgamma(b)-lgamma(c)+lgamma(b+c));
}
double Normalized_PdfParam::operator()(double x)const{
  const unsigned int npar=getNPar();
  double power=pow(x,(*pars[1]))*pow((1-x),(*pars[2]));
  double poly = 1;
  double xx = 1;
  for (unsigned int i = 3; i<npar; i++) {
    xx *= x;
    poly+=(*pars[i])*xx;
  }
  return norm*power*poly;
}
double Normalized_PdfParam::moment(int n)const{
  //Integral of HERAPDF-style function is expressed in terms of Euler beta function:
  //beta(x,y)=int_0^1 t^(x-1)*(1-t)^(y-1) dx=gamma(x)*gamma(y)/gamma(x+y)
  //moment(n)=int_0^1 P[0]*x^(P[1]+n)*(1-x)^P[2]*(1+sum_{i=3}^N{P[i]*x^(i-2)})
  //Let A:=P[0], B:=P[1]+n+1, C:=P[2]+1, then
  //moment=int_0^1 A*x^(B-1)*(1-x)^(C-1)*(1+sum_{i=3}^N{P[i]*x^(i-2)})
  //=A*(beta(B,C)+sum_{i=3}^N{P[i]*beta(B+i-2,C)})
  //beta(B+1,C)=B/(B+C)*beta(B,C)
  //=> beta(B+n,C)=beta(B,C)*product_{k=0}^{n-1}{(B+k)/(B+C+k)}
  //=> beta(B+i-2,C)=beta(B,C)*product_{k=0}^{i-3}{(B+k)/(B+C+k)}
  //=> moment=A*beta(B,C)*(1+sum_{i=3}^N{P[i]*product_{k=0}^{i-3}{(B+k)/(B+C+k)}})=
  //=> moment=A*beta(B,C)*(1+P[3]*B/(B+C)+P[4]*B/(B+C)*(B+1)/(B+C+1)+...)
  //beta(B,C)=exp(lgamma(B)+lgamma(C)-lgamma(B+C))
  const double B0=(*pars[1])+1;
  const double B=B0+n,C=(*pars[2])+1;
  if(B<=0.||C<=0.)return NAN;// integral does not converge
  const size_t N=getNPar();
  double sum=1;
  double prod=1;
  double a=B;
  double b=B+C;
  for(size_t i=3;i<N;++i){
    prod=prod*a/b;
    sum+=(*pars[i])*prod;
    a++;
    b++;
  }
  if(n<0){//beta(x-n,y)/beta(x,y)=product_{k=0}^{n-1}{(x+y+k-n)/(x+k-n)}
    b=B;
    a=B+C;
    const size_t endi=-n;
    for(size_t i=0;i<endi;++i,a+=1,b+=1)sum*=a/b;
  }else if(n>0){//beta(x+n,y)/beta(x,y)=product_{k=0}^{n-1}{(x+k)/(x+y+k)}
    a=B0;
    b=B0+C;
    const size_t endi=n;
    for(size_t i=0;i<endi;++i,a+=1,b+=1)sum*=a/b;
  }
  return (*pars[0])*sum;
}
void Normalized_PdfParam::setMoment(int n,double val){
  BasePdfParam::setMoment(n,val);
  updateNormalization();
}
}
