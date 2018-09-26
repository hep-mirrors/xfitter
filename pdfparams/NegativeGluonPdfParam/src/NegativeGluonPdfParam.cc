 
/*
   @file NegativeGluonPdfParam.cc
   @date 2018-09-26
   @author AddPdfParam.py
   Created by AddPdfParam.py on 2018-09-26
*/

#include "NegativeGluonPdfParam.h"
#include <cmath>

namespace xfitter{
//for dynamic loading
  extern"C" NegativeGluonPdfParam*create(const char*name){
    return new NegativeGluonPdfParam(name);
  }
// Main function to compute PDF
  double NegativeGluonPdfParam::operator()(double x)const{
    //Your code here
    const int npar = getNPar();
    if (npar !=8) {
      return NAN;
    }
    double Pos = pow(x,(*pars[1]))*pow((1-x),(*pars[2])) * ( 1 + x * (*pars[3]) + x*x * (*pars[4])*(*pars[4]));
    double Neg = (*pars[5])*pow(x,(*pars[6]))*pow((1-x),(*pars[7])) ;
    return (*pars[0])*(Pos-Neg);
  }
  double NegativeGluonPdfParam::moment(int n)const{
    /// cut and paste from HERAPDF_PdfParam

    // Positive part:
    const double B=(*pars[1])+(n+1) , C=(*pars[2])+1;
    double sum=1;
    double prod=1;
    double a=B;
    double b=B+C;
    for(int i=3;i<5;++i){
      prod=prod*a/b;
      sum+=(*pars[i])*prod;
      a++;
      b++;
    }
    // Negative part:
    const double Bn=(*pars[6])+(n+1) , Cn=(*pars[7])+1;

    return (*pars[0])*( exp(lgamma(B)+lgamma(C)-lgamma(B+C))*sum - (*pars[5])*exp(lgamma(Bn)+lgamma(Cn)-lgamma(Bn+Cn)) ) ;
  }
  
}
