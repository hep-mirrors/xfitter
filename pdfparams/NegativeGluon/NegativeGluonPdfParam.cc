 
/*
   @file NegativeGluonPdfParam.cc
   @date 2018-09-26
   @author AddPdfParam.py
   Created by AddPdfParam.py on 2018-09-26
*/

#include "NegativeGluonPdfParam.h"
#include "xfitter_cpp_base.h"
#include <cmath>
#include <iostream>

namespace xfitter{
//for dynamic loading
  extern"C" NegativeGluonPdfParam*create(const char*name){
    return new NegativeGluonPdfParam(name);
  }
// Main function to compute PDF
  void NegativeGluonPdfParam::atStart(){
    using namespace std;
    BasePdfParam::atStart();
    const size_t n=getNPar();
    if(n!=8){
      cerr<<"[ERROR] Wrong number of parameters given to parameterisation \""<<_name<<"\", expected 8, got "<<n<<endl;
      hf_errlog(18120700,"F: Wrong number of parameters for a parameterisation, see stderr");
    }
  }
  
  double NegativeGluonPdfParam::operator()(double x) const {
    double Pos = (*pars[0])*pow(x,(*pars[1]))*pow((1-x),(*pars[2])) * ( 1 + x * (*pars[3]) + x*x *(*pars[4]));
    double Neg = (*pars[5])*pow(x,(*pars[6]))*pow((1-x),(*pars[7])) ;
    return Pos - Neg;
  }

  double NegativeGluonPdfParam::moment(int n) const {
    return moment_pos(n) - moment_neg(n);
  }

  double NegativeGluonPdfParam::moment_pos(int n) const {
    /// cut and paste from HERAPDF_PdfParam

    const double B=(*pars[1])+(n+1) , C=(*pars[2])+1;
    if(B<=0.||C<=0.)return NAN;// integral does not converge
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

    return (*pars[0])*( exp(lgamma(B)+lgamma(C)-lgamma(B+C))*sum ) ;
  }

  double NegativeGluonPdfParam::moment_neg(int n) const {
    const double Bn=(*pars[6])+(n+1) , Cn=(*pars[7])+1;
    if(Bn<=0.||Cn<=0.)return NAN;// integral does not converge

    return (*pars[5])*exp(lgamma(Bn)+lgamma(Cn)-lgamma(Bn+Cn)) ;
  }

  void NegativeGluonPdfParam::setMoment(int nMoment,double value) {
    // new_moment = f * pos - neg
    // -> f = (new_moment + neg) / pos
    *pars[0] *= (value + moment_neg(nMoment)) / moment_pos(nMoment);
  }
}
