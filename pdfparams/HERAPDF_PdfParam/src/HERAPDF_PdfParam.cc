/*
   @file HERAPDF_PdfParam.cc
   @date 2018-07-11
   @author  AddPdfParam.py
   Created by  AddPdfParam.py on 2018-07-11
*/

#include "HERAPDF_PdfParam.h"
#include <cmath>

namespace xfitter{
//for dynamic loading
extern"C" HERAPDF_PdfParam*create(const char*name){
  return new HERAPDF_PdfParam(name);
}
// Main function to compute PDF
double HERAPDF_PdfParam::operator()(double x)const{
  const int npar = getNPar();
  if (npar<3) {
    return NAN;
  }
  double power=(*pars[0])*pow(x,(*pars[1]))*pow((1-x),(*pars[2]));
  double poly = 1;
  double xx = 1;
  for (int i = 3; i<npar; i++) {
    xx *= x;
    poly+=(*pars[i])*xx;
  }
  return power*poly;
}
double HERAPDF_PdfParam::moment(int n)const{
  //Integral of HERAPDF-style function is expressed in terms of Euler beta function:
  //beta(x,y)=int_0^1 t^(x-1)*(1-t)^(y-1) dx=gamma(x)*gamma(y)/gamma(x+y)
  //moment(n)=int_0^1 P[0]*x^(P[1]+n)*(1-x)^P[2]*(1+sum_{i=3}^N{P[i]*x^(i-2)})
  //Let A:=P[0], B:=P[1]+n+1, C:=P[2]+1, then
  //moment=int_0^1 A*x^(B-1)*(1-x)^(C-1)*(1+sum_{i=3}^N{P[i]*x^(i-2)})
  //=A*(beta(B,C)+sum_{i=3}^N{P[i]*beta(B+i-2,C)})
  //beta(B+1,C)=B/(B+C)
  //=> beta(B+n,C)=beta(B,C)*product_{k=0}^{n-1}{(B+k)/(B+C+k)}
  //=> beta(B+i-2,C)=beta(B,C)*product_{k=0}^{i-3}{(B+k)/(B+C+k)}
  //=> moment=A*beta(B,C)*(1+sum_{i=3}^N{P[i]*product_{k=0}^{i-3}{(B+k)/(B+C+k)}})=
  //=> moment=A*beta(B,C)*(1+P[3]*B/(B+C)+P[4]*B/(B+C)*(B+1)/(B+C+1)+...)
  //beta(B,C)=exp(lgamma(B)+lgamma(C)-lgamma(B+C))
  using uint=unsigned int;
  const double B=(*pars[1])+(n+1),C=(*pars[2])+1;
  const uint N=getNPar();
  double sum=1;
  double prod=1;
  double a=B;
  double b=B+C;
  for(uint i=3;i<N;++i){
    prod=prod*a/b;
    sum+=(*pars[i])*prod;
    a++;
    b++;
  }
  return (*pars[0])*exp(lgamma(B)+lgamma(C)-lgamma(B+C))*sum;
}
}
