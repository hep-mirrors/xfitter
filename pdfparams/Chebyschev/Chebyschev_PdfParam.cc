#include "Chebyschev_PdfParam.h"
#include "xfitter_cpp_base.h"
#include <cmath>
#include <iostream>


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
  double Chebyschev_PdfParam::operator()(double x) const {
    const unsigned int npar=getNPar();
    double power=(*pars[0])*pow(x,(*pars[1]))*pow((1-x),(*pars[2]));
    
    double poly = 1;
    double xx = 1-2*sqrt(x);
    for (unsigned int i = 3; i<npar; i++) {
      poly+=(*pars[i])*Tn(i-3, xx);
    }
    
    return power*poly;
  }
  
  //Chebyschev polynomial of the first kind
  double Chebyschev_PdfParam::Tn(unsigned int n, double x) const {
    if (n == 0) return 1;
    if (n == 1) return x;
    return 2*x*Tn(n-1, x) - Tn(n-2, x);
  }
  
}

