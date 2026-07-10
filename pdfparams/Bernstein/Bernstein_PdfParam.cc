#include "Bernstein_PdfParam.h"
#include "xfitter_cpp_base.h" //for hf_errlog
#include <cmath>
#include <iostream>
# include "bernstein_polynomial.cpp"


namespace xfitter{

//for dynamic loading
extern"C" Bernstein_PdfParam*create(const char*s){
  return new Bernstein_PdfParam(s);
}


void Bernstein_PdfParam::atStart(){
  using namespace std;
  BasePdfParam::atStart();
  const size_t n=getNPar();
  if(n<3){
    cerr<<"[ERROR] Too few parameters given to parameterisation \""<<_name<<"\", expected at least 3, got "<<n<<endl;
    hf_errlog(18120700,"F: Wrong number of parameters for a parameterisation, see stderr");
  }
}


double Bernstein_PdfParam::operator()(double x)const{
  const unsigned int npar=getNPar();
  double power=(*pars[0])*pow(x,(*pars[1]))*pow((1-x),(*pars[2]));

  double *bernstein_pars=new double[npar-3];
  for (unsigned int i = 3; i<npar; i++) bernstein_pars[i-3] = *pars[i];                                                                                               

  double *power_pars = r8mat_mv_new ( npar-3, npar-3, bernstein_to_power(npar-3), bernstein_pars );

  double poly = 1;
  double xx = 1;
  for (unsigned int i = 3; i<npar; i++) {
    xx *= x;
    poly+=power_pars[i-3]*xx;
  }

  return power*poly;

}



}
  
