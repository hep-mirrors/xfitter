#include "ReactionFractal_DISNC.h"

extern "C" ReactionFractal_DISNC* create() {
  return new ReactionFractal_DISNC();
}

void ReactionFractal_DISNC::compute(TermData*td,valarray<double>&val,map<string,valarray<double> >&err){
  auto&q2=td->getBinColumn("Q2"),
      &x =td->getBinColumn("x"),
      &y =td->getBinColumn("y");

  size_t Npnt=td->getNbins();
  std::valarray<double> f2(Npnt);
  std::valarray<double> fl(Npnt);

  // Get relevant parameters:
  double f_D0  =*td->getParamD("D0");
  double f_D1  =*td->getParamD("D1");
  double f_D2  =*td->getParamD("D2");
  double f_D3  =*td->getParamD("D3");
  double f_Q02 =*td->getParamD("Q02");
  double f_R   =*td->getParamD("R");

  for (size_t i=0; i<Npnt; i++) {

    f2[i] = f_D0 * f_Q02* pow( ( q2[i] / ( q2[i] + f_Q02) ), f_D2-1.0)
    * ( pow(x[i], -f_D2+1.0) ) / ( 1.0 +f_D3 - f_D1 * log (x[i] ))
      * (pow(x[i],-f_D1*log(q2[i]/f_Q02+1.0))
         * pow(q2[i]/f_Q02+1.0,f_D3+1.0)-1.0);
    //    std::cout << f2[i] << " " << pow(q2[i]/f_Q02+1.0,f_D3+1.0) << " " <<  ( pow(x[i], -f_D2+1.0) ) / ( 1.0 +f_D3 - f_D1 * log (x[i] )) << "\n";
  }
  fl = f2*f_R/(1.0+f_R);

  val = f2 - y*y/(1.0+(1.0-y)*(1.0-y))*fl;
}

