 
/*
   @file PolySqrtPdfParam.cc
   @date 2018-08-16
   @author AddPdfParam.py
   Created by AddPdfParam.py on 2018-08-16
*/

#include"PolySqrtPdfParam.h"
#include"xfitter_cpp_base.h"
#include<cmath>
#include<iostream>
using namespace std;
using uint=unsigned int;
namespace xfitter{
//for dynamic loading
extern"C" PolySqrtPdfParam*create(const char*name){
  return new PolySqrtPdfParam(name);
}
void PolySqrtPdfParam::atStart(){
  using namespace std;
  BasePdfParam::atStart();
  const size_t n=getNPar();
  if(n<3){
    cerr<<"[ERROR] Too few parameters given to parameterisation \""<<_name<<"\", expected at least 3, got "<<n<<endl;
    hf_errlog(18120700,"F: Wrong number of parameters for a parameterisation, see stderr");
  }
}
double PolySqrtPdfParam::operator()(double x)const{
  const uint N=getNPar();
  double pol=1;
  double mulx=1;
  double sqrtx=sqrt(x);
  for(uint i=3;i<N;++i){
    mulx*=sqrtx;
    pol+=(*pars[i])*mulx;
  }
  return(*pars[0])*pow(x,(*pars[1]))*pow((1-x),(*pars[2]))*pol;
}
double PolySqrtPdfParam::moment(int n)const{
  //Integral of PolySqrtPdfParam-style function is expressed in terms of Euler beta function:
  //beta(x,y)=int_0^1 t^(x-1)*(1-t)^(y-1) dx=gamma(x)*gamma(y)/gamma(x+y)
  //beta(B,C)=exp(lgamma(B)+lgamma(C)-lgamma(B+C))
  //moment(n)=int_0^1 P[0]*x^(P[1]+n)*(1-x)^P[2]*(1+sum_{i=3}^N{P[i]*x^(i/2-1)})
  //Let A:=P[0], B:=P[1]+n+1, C:=P[2]+1, then
  //moment=int_0^1 A*x^(B-1)*(1-x)^(C-1)*(1+sum_{i=3}^N{P[i]*x^(i/2-1)})
  //=A*(beta(B,C)+sum_{i=3}^N{P[i]*beta(B+i/2-1,C)})=
  //=A*(beta(B,C)+sum(i=3;i<=N;i+=2){P[i]*beta(B+i/2-1,C)}+sum(i=4;i<=N;i+=2){P[i]*beta(B+i/2-1,C)})
  //=A*(beta(B,C)+sum(k=0;3+2k<=N;k++){P[3+2k]*beta(B+k+1/2,C)}+sum(k=1;2+2k<=N;k++){P[2+2k]*beta(B+k,C)})
  //beta(B+1,C)=B/(B+C)
  //=> beta(B+k,C)=beta(B,C)*product_{i=0}^{k-1}{(B+i)/(B+C+i)}
  //=> beta(B+k+1/2,C)=beta(B+1/2,C)*product_{i=0}^{k-1}{(B+1/2+i)/(B+C+1/2+i)}
  //=> moment=A*(beta(B,C)*(1+sum(k=1;2+2k<=N;k++){P[2+2k]*product_{i=0}^{k-1}{(B+i)/(B+C+i)}})+beta(B+1/2,C)*sum(k=0;3+2k<=N;k++){P[3+2k]*product_{i=0}^{k-1}{(B+1/2+i)/(B+C+1/2+i)}})
  using uint=unsigned int;
  const double B=(*pars[1])+(n+1),C=(*pars[2])+1;
  if(B<=0.||C<=0.)return NAN;// integral does not converge
  const uint N=getNPar();
  double sum=1;
  double prod=1;
  double a=B;
  double b=B+C;
  for(uint i=4;i<N;i+=2){
    prod=prod*a/b;
    sum+=(*pars[i])*prod;
    a++;
    b++;
  }
  double lgammaC=lgamma(C);
  double ret=exp(lgamma(B)+lgammaC-lgamma(B+C))*sum;
  sum=0;
  prod=1;
  a=B+0.5;
  b=a+C;
  for(uint i=3;i<N;i+=2){
    sum+=(*pars[i])*prod;
    prod=prod*a/b;
    a++;
    b++;
  }
  ret+=exp(lgamma(B+0.5)+lgammaC-lgamma(B+C+0.5))*sum;
  ret*=(*pars[0]);
  return ret;
}
}
