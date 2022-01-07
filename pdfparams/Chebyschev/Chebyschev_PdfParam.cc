/*
   @file Chebyschev_PdfParam.cc
   @date 2018-07-11
   @author  AddPdfParam.py
   Created by  AddPdfParam.py on 2018-07-11
*/

#include "Chebyschev_PdfParam.h"
#include "xfitter_cpp_base.h"
#include <cmath>
#include <iostream>
#include <boost/math/special_functions/chebyshev.hpp>

namespace xfitter{
//for dynamic loading
extern"C" Chebyschev_PdfParam*create(const char*name){
  return new Chebyschev_PdfParam(name);
}
void Chebyschev_PdfParam::atStart(){
  using namespace std;
  BasePdfParam::atStart();
  const size_t n=getNPar();
  if(n<3){
    cerr<<"[ERROR] Too few parameters given to parameterisation \""<<_name<<"\", expected at least 3, got "<<n<<endl;
    hf_errlog(18120700,"F: Wrong number of parameters for a parameterisation, see stderr");
  }
}
  
// Main function to compute PDF
double Chebyschev_PdfParam::operator()(double x)const{
  const unsigned int npar=getNPar();
  double power=(*pars[0])*pow(x,(*pars[1]))*pow((1-x),(*pars[2]));

  double poly = 1;
  double x = 1-2*sqrt(x);
  for (unsigned int i = 3; i<npar; i++) {
    double Tch = boost::math::chebyshev_t(i, x);
    poly+=(*pars[i])*Tch;
  }
  
  return power*poly;
}
  

}
