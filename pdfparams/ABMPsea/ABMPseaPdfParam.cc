
/*
   @file ABMPseaPdfParam.cc
   @date 2019-02-25
   @author AddPdfParam.py
   Created by AddPdfParam.py on 2019-02-25
*/

#include "ABMPseaPdfParam.h"
#include <cmath>

namespace xfitter{
  //for dynamic loading
  extern"C" ABMPseaPdfParam*create(const char*name){
    return new ABMPseaPdfParam(name);
  }
  // Main function to compute PDF
  double ABMPseaPdfParam::operator()(double x)const{
    const int npar = getNPar();
    if (npar<7) {
      return NAN;
    }
    double power = *pars[1] * (1 + *pars[3] * log(x)) * (1 + *pars[4] * x + *pars[5] * x * x + *pars[6] * x * x * x);
    double val = *pars[0] * pow(1 - x, *pars[2]) * pow(x, power);
    return val;
  }
}
