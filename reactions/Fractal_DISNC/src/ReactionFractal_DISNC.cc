 
/*
   @file ReactionFractal_DISNC.cc
   @date 2017-05-15
   @author  AddReaction.py
   Created by  AddReaction.py on 2017-05-15
*/

#include "ReactionFractal_DISNC.h"

// the class factories
extern "C" ReactionFractal_DISNC* create() {
  return new ReactionFractal_DISNC();
}


// Main function to compute results at an iteration
int ReactionFractal_DISNC::compute(int dataSetID, valarray<double> &val, map<string, valarray<double> > &err)
{
  // Get bin arrays, check that Q2, x and y are present:

  auto *q2p  = GetBinValues(dataSetID,"Q2"), *xp  = GetBinValues(dataSetID,"x"), *yp  = GetBinValues(dataSetID,"y");  
  if (q2p == nullptr || xp == nullptr || yp == nullptr ) {
    std::cout << "\n\nFATAL ERROR: DIS NC requires x,Q2 and y bins to be present !!!" << std::endl;
    std::cout << "CHECK THE DATAFILE !!!" << std::endl;
    return 1;
  }

  auto q2 = *q2p, x = *xp, y = *yp;
  size_t Npnt = q2.size();
  std::valarray<double> f2(Npnt);
  std::valarray<double> fl(Npnt);

  // Get relevant parameters:
  double f_D0  = GetParam("D0");
  double f_Q02 = GetParam("Q02");
  double f_D1  = GetParam("D1");
  double f_D2  = GetParam("D2");
  double f_D3  = GetParam("D3");
  double f_R   = GetParam("R");

  std::valarray<double> part(Npnt);
  for (size_t i=0; i<Npnt; i++) {
    part[i] = pow(x[i],-f_D1*log(q2[i]/f_Q02+1.0));
  }

  f2 = f_D0 * f_Q02* pow( ( q2 / ( q2 + f_Q02) ), f_D2-1.0)
    * ( pow(x, -f_D2+1.0) ) / ( 1+f_D3 - f_D1 * log (x ))
    * part
    * ( pow(q2/f_Q02+1.0,f_D3+1)-1.0);

  fl = f2*f_R/(1.0+f_R);

  val = f2 - y*y/(1.0+(1.0-y)*(1.0-y))*fl;

  std::cout << f_D0 << " " << f_Q02 << " " << f_D1 << " " << f_D2 << " "<< f_D3 << "\n";
  std::cout << q2[0] << " " << f2[0] << " "<< part[0] << std::endl;
  exit(0);

  return 0;
}

