/*
   @file HERAPDF_PdfParam.cc
   @date 2018-07-11
   @author  AddPdfParam.py
   Created by  AddPdfParam.py on 2018-07-11
*/

#include "HERAPDF_PdfParam.h"
#include <cmath>

// Main function to compute PDF
double HERAPDF_PdfParam::compute(double const x, double const* pars)
{
  const int npar = getNPar();
  if (npar<3) {
    return NAN;
  }
  double power = pars[0]*pow(x,pars[1])*pow((1-x),pars[2]);
  double poly = 1;
  double xx = 1;
  for (int i = 3; i<npar; i++) {
    xx *= x;
    poly += pars[i]*xx;
  }
  return power*poly;
}

// Optional function to compute integrals:

// double HERAPDF_PdfParam::moment( double const* pars, int const iMoment)
// {
//   return 0
// }

